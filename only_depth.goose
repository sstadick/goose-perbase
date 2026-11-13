// Difference-counter coverage over indexed BAM; no base unpacking or pileups.
namespace only_depth;
import std;
import dictionary;
import alignment;
import seqair from "github.com/sstadick/goose-seqair@v0.3.1";
import seqair.bam;
import bitvec from "github.com/sstadick/goose-bitvec@v0.1.0";
import opts from "github.com/sstadick/goose-opts@v0.1.0";

struct Options {
    reads:                  opts::Field<const u8[:]> = opts::Field<const u8[:]> { required: true, positional: true, help: "Coordinate-sorted indexed BAM or CRAM 3.x" },
    output:                 opts::Field<const u8[:]> = opts::Field<const u8[:]> { short: "o", help: "Output file; default stdout" },
    ref_fasta:              opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "ref-fasta", short: "r", help: "Optional indexed FASTA; not used for BAM depth" },
    bed_file:               opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "bed-file", short: "b", help: "BED regions of interest" },
    bcf_file:               opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "bcf-file", short: "B", help: "VCF/BCF positions (union with BED)" },
    bgzip:                  opts::Field<bool>        = opts::Field<bool> { short: "Z", help: "BGZF-compress output" },
    compression_threads:    opts::Field<i64>         = opts::Field<i64> { value: 4, long: "compression-threads", short: "T", help: "Output compression workers" },
    compression_level:      opts::Field<i64>         = opts::Field<i64> { value: 2, long: "compression-level", short: "L", help: "Compression level 0..9" },
    threads:                opts::Field<i64>         = opts::Field<i64> { value: hardware_threads(), short: "t", help: "Thread budget; max(threads - 2, 1) workers" },
    chunksize:              opts::Field<i64>         = opts::Field<i64> { value: 1000000, short: "c", help: "Bases per indexed region" },
    channel_size_modifier:  opts::Field<const u8[:]> = opts::Field<const u8[:]> { value: "0.001", long: "channel-size-modifier", short: "C", help: "Approximate output buffering GB per worker" },
    include_flags:          opts::Field<i64>         = opts::Field<i64> { long: "include-flags", short: "f", help: "Require all SAM flag bits" },
    exclude_flags:          opts::Field<i64>         = opts::Field<i64> { long: "exclude-flags", short: "F", help: "Reject any SAM flag bits (default 0)" },
    min_mapq:               opts::Field<i64>         = opts::Field<i64> { long: "min-mapq", short: "q", help: "Minimum mapping quality" },
    mate_fix:               opts::Field<bool>        = opts::Field<bool> { long: "mate-fix", short: "m", help: "Union potentially overlapping same-name mate blocks" },
    fast_mode:              opts::Field<bool>        = opts::Field<bool> { long: "fast-mode", short: "x", help: "Count whole reference spans including skips; conflicts with -m" },
    no_merge:               opts::Field<bool>        = opts::Field<bool> { long: "no-merge", short: "n", help: "Emit individual positions rather than equal-depth runs" },
    bed_format:             opts::Field<bool>        = opts::Field<bool> { long: "bed-format", help: "Headerless five-column BED-like output; use -z for zero-based" },
    skip_merging_intervals: opts::Field<bool>        = opts::Field<bool> { long: "skip-merging-intervals", short: "M", help: "Retain overlapping/duplicate intervals" },
    keep_zeros:             opts::Field<bool>        = opts::Field<bool> { long: "keep-zeros", short: "k", help: "Include uncovered positions with -n; merged runs always include zeros" },
    zero_base:              opts::Field<bool>        = opts::Field<bool> { long: "zero-base", short: "z", help: "Zero-based output positions" },
}

struct Mode {
    enabled:  bool = false,
    fast:     bool = false,
    no_merge: bool = false,
    bed:      bool = false,
}
struct Group {
    name: i64,
    len:  i64,
    next: i64,
}
struct Block {
    group:   i64,
    begin:   i64,
    end:     i64,
    ordinal: i64,
    no_stop: bool,
}

// Coverage changes at each position, relative to the fetched region's start.
// Adding a read increments its start and decrements its end; summing these
// changes gives depth without visiting every covered base in every read.
// With mate fixing, blocks first collect covered intervals by read name so
// overlapping pieces of the same pair contribute only once.
root struct Counter {
    delta:  i64[>..]             = [],
    blocks: Block[>..]           = [],
    names:  u8[>..]              = [],
    groups: Group[>..]           = [],
    ids:    dictionary<u64, i64> = dictionary<u64, i64> {},
}

fn group(counter: Counter&, name: const u8[:]) -> i64 {
    var hash: u64 = 14695981039346656037;
    for c in name {
        hash = (hash ^ c) * 1099511628211;
    }
    let head = counter.ids.get(hash);
    let next = if head { head } else { -1 };
    var id = next;
    while id >= 0 {
        let g = counter.groups[id];
        if name == counter.names[g.name..g.name + g.len] {
            return id;
        }
        id = g.next;
    }
    id = counter.groups.len;
    counter.groups.push(Group { name: counter.names.len, len: name.len, next: next });
    counter.names.append(name);
    counter.ids.insert(hash, id);
    id
}

fn add(counter: Counter&, changed: bitvec::BitVec&, begin: i64, end: i64, no_stop: bool) {
    counter.delta[begin] += 1;
    changed.set(begin);
    if !no_stop {
        counter.delta[end] -= 1;
        changed.set(end);
    }
}

fn add_block(counter: Counter&, changed: bitvec::BitVec&, begin: i64, end: i64, start: i64, stop: i64, group: i64) {
    if begin >= stop || end <= start {
        return;
    }
    let a = max(begin, start) - start;
    let b = min(end - start, stop - start - 1);
    if group < 0 {
        counter.add(changed, a, b, end >= stop);
    } else {
        counter.blocks.push(Block {
            group:   group,
            begin:   a,
            end:     b,
            ordinal: counter.blocks.len,
            no_stop: end >= stop,
        });
    }
}

fn load<R, I, C>(counter: Counter&, changed: bitvec::BitVec&, reader: R&, index: I&, cursor: C&, source: alignment::Reference&,
        tid: i32, start: i64, stop: i64,
        include: u16, exclude: u16, mapq: u8, mate_fix: bool, fast: bool) -> bam::Status {
    // The bit vector marks counter positions written by the previous region.
    // Reset only those positions, before resizing for a possibly shorter region.
    changed.ones() { pos =>
        counter.delta[pos] = 0;
    };
    counter.delta.resize(stop - start, 0);
    changed.reset(stop - start);
    counter.blocks.clear();
    counter.names.clear();
    counter.groups.clear();
    counter.ids.clear();
    let query = alignment::fetch(reader, index, cursor, tid, start, stop);
    if query != bam::Status.Record {
        return query;
    }
    loop {
        let status = alignment::read(reader, cursor, source);
        if status == bam::Status.End {
            break;
        }
        if status != bam::Status.Record {
            return status;
        }
        let record = reader.record();
        let f = record.fields;
        if f.flags & include != include || f.flags & exclude != 0 || f.mapq < mapq {
            continue;
        }
        // HTSlib gives zero-reference-span alignments a one-base fetch span.
        if fast {
            counter.add_block(changed, f.pos, f.end, start, stop, -1);
            continue;
        }
        let eligible = mate_fix && f.flags & 1 != 0 && f.tid == f.mate_tid &&
            (abs(f.template_len as i64) / 2 < f.end - f.pos || abs(f.mate_pos - f.pos) < f.end - f.pos);
        let id = if eligible { counter.group(record.name()) } else { -1 };
        var pos = f.pos;
        for i in f.cigar_len {
            let word = record.cigar(i);
            let op = word & 15;
            let n = (word >> 4) as i64;
            if op == 0 || op == 2 || op == 7 || op == 8 {
                counter.add_block(changed, pos, pos + n, start, stop, id);
                pos += n;
            } else if op == 3 {
                pos += n;
            }
        }
    }
    sort(counter.blocks) { a, b => a.group < b.group || (a.group == b.group &&
        (a.begin < b.begin || (a.begin == b.begin && (a.end < b.end || (a.end == b.end && a.ordinal < b.ordinal))))) };
    var i: i64 = 0;
    while i < counter.blocks.len {
        var current = counter.blocks[i];
        i += 1;
        while i < counter.blocks.len && counter.blocks[i].group == current.group && counter.blocks[i].begin <= current.end {
            current.end = max(current.end, counter.blocks[i].end);
            i += 1;
            // Perbase/rust-lapper retains the first interval's value, even when
            // a later clipped interval extends to the region's final base.
        }
        counter.add(changed, current.begin, current.end, current.no_stop);
    }
    bam::Status.Record
}

fn row(text: u8[>..]&, name: const u8[:], begin: i64, end: i64, depth: i64, bed: bool) {
    var quoted = false;
    for c in name {
        if c == '\t' || c == '\r' || c == '\n' || c == '"' {
            quoted = true;
        }
    }
    if quoted {
        text.push('"');
        for c in name {
            if c == '"' {
                text.push('"');
            }
            text.push(c);
        }
        text.push('"');
    } else {
        text.append(name);
    }
    format(text, "\t", begin, "\t", end, if bed { "\t\t" } else { "\t" }, depth, "\n");
}

// Format coverage runs, invoking F after each complete row so the caller can
// flush text without retaining the whole region. Splitting the output buffer
// must not split or otherwise change the coverage runs themselves.
fn write<F: fn()>(counter: const Counter&, changed: const bitvec::BitVec&, text: u8[>..]&, name: const u8[:], start: i64, base: i64, mode: Mode, zeros: bool) {
    if counter.delta.len == 0 {
        return;
    }
    var depth: i64 = 0;
    var begin: i64 = 0;
    fn emit(end: i64) {
        if mode.no_merge {
            if depth != 0 || zeros {
                for pos in begin..end {
                    row(text, name, start + pos + base, start + pos + base + 1, depth, mode.bed);
                    F();
                }
            }
        } else {
            row(text, name, start + begin + base, start + end + base, depth, mode.bed);
            F();
        }
    }
    changed.ones() { pos =>
        let next = depth + counter.delta[pos];
        if next != depth {
            if pos > begin {
                emit(pos);
            }
            begin = pos;
            depth = next;
        }
    };
    emit(counter.delta.len);
}
