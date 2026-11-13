import std;
import os;
import perbase;
import output;
import result_stream;
import regions;
import only_depth;
import merge_adjacent;
import alignment;
import bitvec from "github.com/sstadick/goose-bitvec@v0.1.0";
import opts from "github.com/sstadick/goose-opts@v0.1.0";
import seqair from "github.com/sstadick/goose-seqair@v0.4.0";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import seqair.cram;
import seqair.cram_index;
import bio from "github.com/sstadick/goose-bio@v0.4.1";
import bio.faidx;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;

struct Options {
    reads:                    opts::Field<const u8[:]> = opts::Field<const u8[:]> { required: true, positional: true, help: "Coordinate-sorted indexed BAM or CRAM 3.x" },
    output:                   opts::Field<const u8[:]> = opts::Field<const u8[:]> { short: "o", help: "Output file; default stdout" },
    ref_fasta:                opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "ref-fasta", short: "r", help: "FASTA with .fai; BGZF also requires .gzi" },
    bed_file:                 opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "bed-file", short: "b", help: "BED regions of interest" },
    bcf_file:                 opts::Field<const u8[:]> = opts::Field<const u8[:]> { long: "bcf-file", short: "B", help: "VCF/BCF positions (union with BED)" },
    bgzip:                    opts::Field<bool>        = opts::Field<bool> { short: "Z", help: "BGZF-compress output" },
    compression_threads:      opts::Field<i64>         = opts::Field<i64> { value: 4, long: "compression-threads", short: "T", help: "Output compression workers" },
    compression_level:        opts::Field<i64>         = opts::Field<i64> { value: 2, long: "compression-level", short: "L", help: "Compression level 0..9" },
    threads:                  opts::Field<i64>         = opts::Field<i64> { value: hardware_threads(), short: "t", help: "Thread budget; max(threads - 2, 1) pileup workers" },
    chunksize:                opts::Field<i64>         = opts::Field<i64> { value: 1000000, short: "c", help: "Bases per indexed region" },
    channel_size_modifier:    opts::Field<const u8[:]> = opts::Field<const u8[:]> { value: "0.15", long: "channel-size-modifier", short: "C", help: "Approximate output buffering GB per worker" },
    include_flags:            opts::Field<i64>         = opts::Field<i64> { long: "include-flags", short: "f", help: "Require all SAM flag bits" },
    exclude_flags:            opts::Field<i64>         = opts::Field<i64> { long: "exclude-flags", short: "F", help: "Reject any SAM flag bits (default 0)" },
    mate_fix:                 opts::Field<bool>        = opts::Field<bool> { long: "mate-fix", short: "m", help: "Resolve overlapping query-name groups" },
    mate_resolution_strategy: opts::Field<const u8[:]> = opts::Field<const u8[:]> { value: "original", long: "mate-resolution-strategy", help: "Original, N, IUPAC, or BaseQual/MapQual priority strategies" },
    skip_merging_intervals:   opts::Field<bool>        = opts::Field<bool> { long: "skip-merging-intervals", short: "M", help: "Retain overlapping/duplicate intervals" },
    min_mapq:                 opts::Field<i64>         = opts::Field<i64> { long: "min-mapq", short: "q", help: "Minimum mapping quality" },
    min_base_quality_score:   opts::Field<i64>         = opts::Field<i64> { long: "min-base-quality-score", short: "Q", help: "Count lower-quality bases as N" },
    keep_zeros:               opts::Field<bool>        = opts::Field<bool> { long: "keep-zeros", short: "k", help: "Include uncovered positions" },
    zero_base:                opts::Field<bool>        = opts::Field<bool> { long: "zero-base", short: "z", help: "Zero-based output positions" },
    ref_cache_size:           opts::Field<i64>         = opts::Field<i64> { value: 10, long: "ref-cache-size", help: "Maximum cached reference windows per worker" },
    max_depth:                opts::Field<i64>         = opts::Field<i64> { value: 100000, long: "max-depth", short: "D", help: "HTSlib-compatible read admission limit" },
    max_region_bytes:         opts::Field<i64>         = opts::Field<i64> { value: 4294967295, long: "max-region-bytes", help: "Retained payload bytes per region, 1..4294967295; not total RAM" },
}

// Each worker receives its own copy of these options and opens its own input
// files. No BAM/CRAM reader or borrowed record crosses a thread boundary.
struct Config {
    reads:            u8[],
    index:            u8[],
    reference:        u8[],
    filters:          perbase::Filters,
    mate_mode:        i64,
    keep_zeros:       bool,
    coordinate_base:  i64,
    max_depth:        i64,
    cache_size:       i64,
    max_region_bytes: i64,
    depth:            only_depth::Mode,
}
struct Job {
    id:    i64,
    tid:   i32,
    begin: i64,
    end:   i64,
}
struct CacheEntry {
    tid:    i32,
    begin:  i64,
    end:    i64,
    offset: i64,
}
root struct ReferenceCache {
    entries: CacheEntry[>..] = [],
    bytes:   u8[>..]         = [],
}

fn fail(message: const u8[:]) {
    write_stderr(str("goose-perbase: ", message, "\n"));
    exit(1);
}

fn check(ok: bool, message: const u8[:]) {
    if !ok {
        fail(message);
    }
}

fn header(text: u8[>..]&, reference: bool) {
    text.append("REF\tPOS\t");
    if reference {
        text.append("REF_BASE\t");
    }
    text.append("DEPTH\tA\tC\tG\tT\tN\tR\tY\tS\tW\tK\tM\tINS\tDEL\tREF_SKIP\tFAIL\tCOUNT_OF_MATE_RESOLUTIONS\tNEAR_MAX_DEPTH\n");
}

fn fetch(cache: ReferenceCache&, reference: faidx::Reader&, name: const u8[:], job: Job, capacity: i64, bases: u8[>..]&) {
    bases.clear();
    var replace = cache.entries.len >= capacity;
    for entry in cache.entries {
        if entry.tid == job.tid {
            if entry.begin <= job.begin && entry.end >= job.end {
                let start = entry.offset + job.begin - entry.begin;
                bases.append(cache.bytes[start..start + job.end - job.begin]);
                return;
            }
            replace = true;
        }
    }
    check(reference.fetch(name, job.begin, job.end, bases) == faidx::Status.Ok, "reference range does not match BAM");
    if capacity == 0 {
        return;
    }
    if replace {
        cache.entries.clear();
        cache.bytes.clear();
    }
    cache.entries.push(CacheEntry { tid: job.tid, begin: job.begin, end: job.end, offset: cache.bytes.len });
    cache.bytes.append(bases);
}

thread_fn worker(config: Config, stream: u64, window: i64) {
    let bam_reader = bam::Reader { .. };
    let index = bam_index::Index { .. };
    let reference = faidx::Reader { .. };
    let cram_reader = cram::Reader { .. };
    let crai = cram_index::Index { .. };
    let cursor = cram_index::Cursor { .. };
    defer bam_reader.close();
    defer cram_reader.close();
    defer reference.close();
    let reader = alignment::Header { .. };
    let bam_cursor = alignment::BamCursor { .. };
    let source = alignment::Reference { reader: reference };
    check(reader.open(bam_reader, index, cram_reader, crai, reference, config.reads, config.index, config.reference) == bam::Status.Ok, "cannot open indexed alignment or reference");
    let has_reference = config.reference.len > 0;
    // Counting never reads auxiliary tags or Seqair's primary-pair links.
    // Mate resolution below groups all same-name reads using its own rules.
    let store = pileup::Store {
        max_bytes:        config.max_region_bytes,
        retain_auxiliary: false,
        auto_link_mates:  false,
        ..
    };
    let mates = perbase::MateCounter { .. };
    let cache = ReferenceCache { .. };
    let depth = only_depth::Counter { .. };
    let changed = bitvec::BitVec { .. };
    var text:     u8[>..]            = [];
    var bases:    u8[>..]            = [];
    var admitted: bool[>..]          = [];
    var selected: pileup::Entry[>..] = [];
    loop {
        let job = qget<Job>();
        if job.id < 0 {
            break;
        }
        let slot = job.id % window;
        text.clear();
        if config.depth.enabled {
            let status = if reader.is_cram {
                depth.load(changed, cram_reader, crai, cursor, source, job.tid, job.begin, job.end, config.filters.include_flags,
                    config.filters.exclude_flags, config.filters.min_mapq, config.mate_mode > 0, config.depth.fast)
            } else {
                depth.load(changed, bam_reader, index, bam_cursor, source, job.tid, job.begin, job.end, config.filters.include_flags,
                    config.filters.exclude_flags, config.filters.min_mapq, config.mate_mode > 0, config.depth.fast)
            };
            check(status == bam::Status.Ok, "cannot read indexed alignment region (check reference for CRAM)");
            let ref = reader.references[job.tid as i64];
            let name = reader.names[ref.name..ref.name + ref.name_len];
            depth.write(changed, text, name, job.begin, config.coordinate_base, config.depth, config.keep_zeros) {
                result_stream::flush(stream, slot, text);
            };
            result_stream::flush(stream, slot, text, true);
            continue;
        }
        let status = if reader.is_cram { alignment::load(store, cram_reader, crai, cursor, source, job.tid, job.begin, job.end) }
            else { store.load_indexed(bam_reader, index, job.tid, job.begin, job.end,
                pileup::Filter { exclude_flags: 0, min_mapq: 0, retain_zero_span: true }) };
        check(status == bam::Status.Ok, "cannot load alignment region (reference, invalid data, sorting, or size limit; --max-region-bytes controls retained payload)");
        let restricted = perbase::admit(store, min(config.max_depth, 2147483647), admitted);
        if config.mate_mode > 0 {
            mates.prepare(store);
        }
        let ref = reader.references[job.tid as i64];
        let name = reader.names[ref.name..ref.name + ref.name_len];
        if has_reference && (store.records.len > 0 || config.keep_zeros) {
            cache.fetch(reference, name, job, config.cache_size, bases);
        }
        var next = job.begin;
        store.columns() { position, raw =>
            if restricted {
                selected.clear();
                for entry in raw {
                    if admitted[entry.record_index] {
                        selected.push(entry);
                    }
                }
            }
            let entries = if restricted { selected[..] } else { raw };
            if entries.len > 0 {
                if config.keep_zeros {
                    while next < position {
                        perbase::append_row(text, name, next + config.coordinate_base, if has_reference { bases[next - job.begin] } else { 0 }, perbase::Counts { .. }, false);
                        result_stream::flush(stream, slot, text);
                        next += 1;
                    }
                }
                let counts = if config.mate_mode > 0 { mates.count(store, entries, config.filters, config.mate_mode) } else { perbase::count(store, entries, config.filters) };
                perbase::append_row(text, name, position + config.coordinate_base, if has_reference { bases[position - job.begin] } else { 0 }, counts,
                    entries.len > config.max_depth - config.max_depth / 100);
                result_stream::flush(stream, slot, text);
                next = position + 1;
            }
        };
        if config.keep_zeros {
            while next < job.end {
                perbase::append_row(text, name, next + config.coordinate_base, if has_reference { bases[next - job.begin] } else { 0 }, perbase::Counts { .. }, false);
                result_stream::flush(stream, slot, text);
                next += 1;
            }
        }
        result_stream::flush(stream, slot, text, true);
    }
}

fn main() {
    let argv = args();
    if argv.len < 2 || argv[1] == "--help" || argv[1] == "-h" {
        write_stdout("Usage: goose-perbase <base-depth|only-depth|merge-adjacent> [OPTIONS] [INPUT]\n");
        return;
    }
    if argv[1] == "merge-adjacent" {
        let merge_options, help, error = opts::parse<merge_adjacent::Options>(argv[2..]);
        if error.len != 0 {
            fail(error);
        }
        if help {
            write_stdout(opts::help<merge_adjacent::Options>("goose-perbase merge-adjacent"));
            return;
        }
        merge_adjacent::run(merge_options);
        return;
    }
    if argv[1] == "only-depth" {
        run_only(argv[2..]);
        return;
    }
    if argv[1] != "base-depth" {
        fail("expected base-depth, only-depth or merge-adjacent");
    }
    run_base(argv[2..]);
}

fn run_base<S>(argv: S[:]) {
    if argv.len == 0 {
        write_stdout(opts::help<Options>("goose-perbase base-depth"));
        return;
    }
    let options, help, error = opts::parse<Options>(argv);
    if error.len != 0 {
        fail(error);
    }
    if help {
        write_stdout(opts::help<Options>("goose-perbase base-depth"));
        return;
    }
    run(options);
}

fn run_only<S>(argv: S[:]) {
    if argv.len == 0 {
        write_stdout(opts::help<only_depth::Options>("goose-perbase only-depth"));
        return;
    }
    let d, help, error = opts::parse<only_depth::Options>(argv);
    if error.len != 0 {
        fail(error);
    }
    if help {
        write_stdout(opts::help<only_depth::Options>("goose-perbase only-depth"));
        return;
    }
    check(!d.mate_fix.value || !d.fast_mode.value, "mate-fix conflicts with fast-mode");
    let common = Options {
        reads:                  copy(d.reads),
        output:                 copy(d.output),
        ref_fasta:              copy(d.ref_fasta),
        bed_file:               copy(d.bed_file),
        bcf_file:               copy(d.bcf_file),
        bgzip:                  copy(d.bgzip),
        compression_threads:    copy(d.compression_threads),
        compression_level:      copy(d.compression_level),
        threads:                copy(d.threads),
        chunksize:              copy(d.chunksize),
        channel_size_modifier:  copy(d.channel_size_modifier),
        include_flags:          copy(d.include_flags),
        exclude_flags:          copy(d.exclude_flags),
        mate_fix:               copy(d.mate_fix),
        skip_merging_intervals: copy(d.skip_merging_intervals),
        min_mapq:               copy(d.min_mapq),
        keep_zeros:             copy(d.keep_zeros),
        zero_base:              copy(d.zero_base),
        ..
    };
    run(common, only_depth::Mode {
        enabled:  true,
        fast:     d.fast_mode.value,
        no_merge: d.no_merge.value,
        bed:      d.bed_format.value,
    });
}

fn run(options: const Options&, depth: only_depth::Mode = only_depth::Mode { .. }) {
    check(options.threads.value > 0 && options.compression_threads.value > 0, "thread counts must be positive");
    check(options.chunksize.value > 0 && options.chunksize.value <= 4294967295, "invalid chunksize");
    check(options.max_depth.value >= 0 && options.max_depth.value <= 4294967295, "invalid max-depth");
    check(options.max_region_bytes.value > 0 && options.max_region_bytes.value <= 4294967295, "max-region-bytes must be 1..4294967295");
    check(options.ref_cache_size.value >= 0 && options.ref_cache_size.value <= 1048576, "invalid ref-cache-size");
    check(options.compression_level.value >= 0 && options.compression_level.value <= 9, "invalid compression-level");
    let modifier, modifier_ok = parse_flt(options.channel_size_modifier.value);
    check(modifier_ok && modifier > 0 && modifier <= 1000000, "invalid channel-size-modifier");
    check(options.include_flags.value >= 0 && options.include_flags.value <= 65535 && options.exclude_flags.value >= 0 && options.exclude_flags.value <= 65535 &&
        options.min_mapq.value >= 0 && options.min_mapq.value <= 255 && options.min_base_quality_score.value >= 0 && options.min_base_quality_score.value <= 255, "filter out of range");
    for path in [options.reads.value, options.ref_fasta.value, options.bed_file.value, options.bcf_file.value] {
        check(options.output.value.len == 0 || options.output.value == "-" || options.output.value != path, "output would overwrite input");
    }
    let mode = perbase::strategy(options.mate_resolution_strategy.value);
    check(mode > 0, "unknown mate-resolution-strategy");
    let bam_reader = bam::Reader { .. };
    let index = bam_index::Index { .. };
    let reference = faidx::Reader { .. };
    let cram_reader = cram::Reader { .. };
    let crai = cram_index::Index { .. };
    let cursor = cram_index::Cursor { .. };
    defer bam_reader.close();
    defer cram_reader.close();
    defer reference.close();
    let reader = alignment::Header { .. };
    check(reader.open(bam_reader, index, cram_reader, crai, reference, options.reads.value, "", options.ref_fasta.value) == bam::Status.Ok, "cannot open indexed alignment or reference");
    let workers = max(1, min(options.threads.value, hardware_threads()) - 2);
    let config = Config {
        reads:     options.reads.value,
        index:     reader.index_path[..],
        reference: options.ref_fasta.value,
        filters: perbase::Filters {
            include_flags:    options.include_flags.value as u16,
            exclude_flags:    options.exclude_flags.value as u16,
            min_mapq:         options.min_mapq.value as u8,
            min_base_quality: options.min_base_quality_score.value as u8,
        },
        mate_mode:        if options.mate_fix.value { mode } else { 0 },
        keep_zeros:       options.keep_zeros.value,
        coordinate_base:  if options.zero_base.value { 0 } else { 1 },
        max_depth:        options.max_depth.value,
        cache_size:       options.ref_cache_size.value,
        max_region_bytes: options.max_region_bytes.value,
        depth:            depth,
    };
    var intervals: regions::Interval[>..] = [];
    if options.bed_file.value.len > 0 {
        regions::read(options.bed_file.value, false, reader, intervals);
    }
    if options.bcf_file.value.len > 0 {
        regions::read(options.bcf_file.value, true, reader, intervals);
    }
    let restricted = options.bed_file.value.len > 0 || options.bcf_file.value.len > 0;
    if restricted {
        regions::merge(intervals, !options.skip_merging_intervals.value);
    } else {
        for ref, tid in reader.references {
            var begin: i64 = 0;
            while begin < ref.length {
                let end = min(ref.length, begin + options.chunksize.value);
                intervals.push(regions::Interval { tid: tid as i32, begin: begin, end: end, .. });
                begin = end;
            }
        }
    }
    var jobs:         Job[>..] = [];
    var longest_name: i64      = 0;
    let width = min(4294967295, options.chunksize.value * workers);
    var first: i64 = 0;
    for ref, tid in reader.references {
        let name = reader.names[ref.name..ref.name + ref.name_len];
        longest_name = max(longest_name, name.len);
        var begin: i64 = 0;
        while begin < ref.length {
            let end = min(ref.length, begin + width);
            while first < intervals.len && (intervals[first].tid < tid || (intervals[first].tid == tid && intervals[first].end <= begin)) {
                first += 1;
            }
            var i = first;
            while i < intervals.len && intervals[i].tid == tid && intervals[i].begin < end {
                if intervals[i].end > begin {
                    jobs.push(Job {
                        id:    jobs.len,
                        tid:   tid as i32,
                        begin: max(begin, intervals[i].begin),
                        end:   min(end, intervals[i].end),
                    });
                }
                i += 1;
            }
            begin = end;
        }
    }
    // Keep the existing -C concurrency policy and the exact genomic jobs.
    // Changing region boundaries would change maximum-depth read admission.
    let estimate = if depth.enabled { 65536 } else { max(1, options.chunksize.value * (longest_name + 64)) };
    let budget = modifier * 1000000000.0 * workers as f64;
    // Round down to whole jobs. The value is bounded by the worker count;
    // discarding its fractional part is intentional, including in checked builds.
    let window = min(workers, max(1, (min(budget / estimate as f64, workers as f64)) as! i64));
    // Base-depth needs room to count and format ahead of the ordered writer.
    // A single-piece mailbox makes later workers stop almost immediately.
    // Allow up to one estimated region per mailbox within the -C budget,
    // with a 64 KiB minimum so even a tiny budget can make progress.
    let capacity = max(65536, (min(estimate as f64, budget / window as f64)) as! i64);
    qlimit<Job>(workers);
    let stream = result_stream::open(window, capacity);
    check(stream != 0, "cannot allocate result mailboxes");
    defer result_stream::close(stream);
    var ids: i64[>..] = [];
    for i in min(window, jobs.len) {
        ids.push(thread_spawn(worker, copy(config), stream, window));
    }
    let sink = output::Sink { .. };
    sink.open(options.output.value, options.bgzip.value, min(options.compression_threads.value, hardware_threads()), options.compression_level.value as i32);
    var received: u8[65536] = [0; 65536];
    var started = false;
    var submitted = min(window, jobs.len);
    for i in submitted {
        check(qput(jobs[i]), "job queue full");
    }
    for at in jobs.len {
        loop {
            let size = result_stream::receive(stream, at % window, received);
            if size == 0 {
                break;
            }
            if !started {
                var text: u8[>..] = [];
                if depth.enabled {
                    if !depth.bed {
                        text.append("REF\tPOS\tEND\tDEPTH\n");
                    }
                } else {
                    header(text, config.reference.len > 0);
                }
                sink.write(text);
                started = true;
            }
            sink.write(received[..size]);
        }
        // Refill as soon as this region is written, rather than waiting for
        // the rest of a batch. Only this now-empty mailbox may be reused.
        if submitted < jobs.len {
            check(qput(jobs[submitted]), "job queue full");
            submitted += 1;
        }
    }
    for id in ids {
        check(qput(Job { id: -1, tid: -1, begin: 0, end: 0 }), "job queue full");
    }
    for id in ids {
        check(thread_wait(id), "worker failed");
    }
    sink.close();
}
