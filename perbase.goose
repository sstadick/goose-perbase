// Count bases, deletions, and other alignment events at one reference position.
// Each input entry describes how one read overlaps that position.
namespace perbase;
import std;
import seqair from "github.com/sstadick/goose-seqair@v0.4.0";
import seqair.bam;
import seqair.pileup;

// Look up a mate-resolution strategy by its case-insensitive CLI name.
// Return a positive strategy ID, or -1 for an unknown name. The CLI uses zero
// separately to mean that overlapping reads should not be combined.
fn strategy(name: const u8[:]) -> i64 {
    let names = ["original", "basequalmapqualfirstinpair", "basequalmapqualiupac",
        "basequalmapqualn", "mapqualbasequalfirstinpair", "mapqualbasequaliupac",
        "mapqualbasequaln", "iupac", "n"];
    for candidate, i in names {
        if candidate.len != name.len {
            continue;
        }
        var equal = true;
        for c, j in name {
            if (if c >= 'A' && c <= 'Z' { c + 32 } else { c }) != candidate[j] {
                equal = false;
            }
        }
        if equal {
            return i + 1;
        }
    }
    -1
}

// Rules for counting reads at a reference position. include_flags requires
// every specified SAM flag bit; exclude_flags rejects any specified bit.
// min_mapq is the minimum read mapping quality. Defaults accept all flag and
// mapping-quality values, matching Rust Perbase.
//
// Keep rejected reads when loading the region: counting records them in FAIL
// rather than DEPTH. min_base_quality is different: a passing read with a base
// below this quality contributes N and still contributes to DEPTH.
//
// ```goose
// let options = perbase::Filters { .. };
// assert(options.exclude_flags == 0 && options.min_base_quality == 0);
// ```
struct Filters {
    include_flags:    u16 = 0,
    exclude_flags:    u16 = 0,
    min_mapq:         u8  = 0,
    min_base_quality: u8  = 0,
}

// Counts for one reference position, in output-column order:
// DEPTH, A, C, G, T, N, R, Y, S, W, K, M, INS, DEL, REF_SKIP, FAIL,
// COUNT_OF_MATE_RESOLUTIONS. values[0] is DEPTH and values[16] is the number
// of passing read-name groups where multiple alignments were resolved.
//
// CIGAR deletions contribute to DEPTH and DEL. Reference skips, such as introns,
// contribute to REF_SKIP but not DEPTH. INS counts insertion events anchored
// after this position, not the number of inserted bases. Ambiguous bases in
// input reads count as N; mate resolution can produce the R/Y/S/W/K/M counts.
//
// ```goose
// let row = perbase::Counts { .. };
// assert(row.values[0] == 0 && row.values[16] == 0);
// ```
struct Counts {
    values: i64[17] = [0; 17],
}

// Count one reference position without combining overlapping mates.
// Pass the entries supplied by store.columns() and the same store that
// produced them. Each entry identifies one read and its base, deletion, or
// reference skip at this position. The returned Counts owns only numbers
// and remains valid after the store is reused.
//
// Missing BAM base qualities are represented by 255. That value participates
// in the ordinary numeric quality comparison, matching Rust Perbase.
//
// ```goose
// let store = pileup::Store { .. };
// let row = perbase::count(store, [], perbase::Filters { .. });
// assert(row.values[0] == 0);
// ```
fn count(
    store:   const pileup::Store&,
    entries: const pileup::Entry[:],
    filters: Filters,
) -> Counts {
    var row = Counts { .. };
    for &entry in entries {
        add(row, entry, filters, passes(store, entry, filters), 0, false);
    }
    row
}

fn passes(store: const pileup::Store&, entry: pileup::Entry, filters: Filters) -> bool {
    let f = store.records[entry.record_index];
    (f.flags & filters.include_flags) == filters.include_flags && f.flags & filters.exclude_flags == 0 &&
        f.mapq >= filters.min_mapq
}

checkfree inline fn add(
    row: Counts&,
    entry: pileup::Entry,
    filters: Filters,
    passed: bool,
    recommended: u8,
    resolved: bool,
) {
    if !passed {
        row.values[15] += 1;
        return;
    }
    if resolved {
        row.values[16] += 1;
    }
    if entry.op == pileup::Op.RefSkip {
        row.values[14] += 1;
        return;
    }
    row.values[0] += 1;
    if entry.op == pileup::Op.Deletion {
        row.values[13] += 1;
        return;
    }
    var base = if recommended != 0 { recommended } else { entry.base };
    if entry.quality < filters.min_base_quality {
        base = 'N';
    }
    if base == 'A' {
        row.values[1] += 1;
    } else if base == 'C' {
        row.values[2] += 1;
    } else if base == 'G' {
        row.values[3] += 1;
    } else if base == 'T' || base == 'U' {
        row.values[4] += 1;
    } else if recommended != 0 {
        var found = false;
        let ambiguity: u8[6] = ['R', 'Y', 'S', 'W', 'K', 'M'];
        const for <I> in [0, 1, 2, 3, 4, 5] {
            if ambiguity[I] == base {
                row.values[6 + I] += 1;
                found = true;
            }
        }
        if !found {
            row.values[5] += 1;
        }
    } else {
        row.values[5] += 1;
    }
    if entry.indel > 0 {
        row.values[12] += 1;
    }
}

fn canonical(base: u8) -> u8 {
    if base == 'A' || base == 'C' || base == 'G' || base == 'T' {
        base
    } else if base == 'U' {
        'T'
    } else {
        'N'
    }
}

inline fn choose(
    store: const pileup::Store&,
    a: pileup::Entry,
    b: pileup::Entry,
    ap: bool,
    bp: bool,
    mode: i64,
) -> bool, u8 {
    if ap != bp {
        return ap, 0;
    }
    if !ap {
        return true, 0;
    }
    let af = store.records[a.record_index];
    let bf = store.records[b.record_index];
    let fallback = mode == 1 || a.qpos < 0 || b.qpos < 0;
    if fallback {
        if af.mapq != bf.mapq {
            return af.mapq > bf.mapq, 0;
        }
    } else if mode < 8 {
        let aq = if mode <= 4 { a.quality } else { af.mapq };
        let bq = if mode <= 4 { b.quality } else { bf.mapq };
        if aq != bq {
            return aq > bq, 0;
        }
        let ar = if mode <= 4 { af.mapq } else { a.quality };
        let br = if mode <= 4 { bf.mapq } else { b.quality };
        if ar != br {
            return ar > br, 0;
        }
    }
    if fallback || mode == 2 || mode == 5 {
        return af.flags & 64 != 0 || bf.flags & 64 == 0, 0;
    }
    let ab = canonical(a.base);
    let bb = canonical(b.base);
    if ab == bb {
        return true, ab;
    }
    if mode == 4 || mode == 7 || mode == 9 {
        return true, 'N';
    }
    let lo = min(ab, bb);
    let hi = max(ab, bb);
    var base: u8 = 'N';
    if lo == 'A' && hi == 'G' {
        base = 'R';
    } else if lo == 'C' && hi == 'T' {
        base = 'Y';
    } else if lo == 'C' && hi == 'G' {
        base = 'S';
    } else if lo == 'A' && hi == 'T' {
        base = 'W';
    } else if lo == 'G' && hi == 'T' {
        base = 'K';
    } else if lo == 'A' && hi == 'C' {
        base = 'M';
    }
    return true, base;
}

struct ReadInfo {
    group: i64,
    passed: bool = true,
}

// Reusable storage for resolving reads with the same name at each position.
// prepare() groups the region's records by read name, including secondary
// alignments and unpaired reads. This matches Perbase's grouping rule rather
// than restricting groups to primary mate pairs. Coordinate-sorted records
// whose spans cannot overlap are counted directly at each position.
root struct MateCounter {
    groups:  ReadInfo[>..] = [],
    heads:   i64[>..] = [],
    links:   i64[>..] = [],
    touched: i64[>..] = [],
    members: i64[>..] = [],
    group_ends: i64[>..] = [],
    overlaps: bool[>..] = [],
    table: i64[>..] = [],
    filters_ready: bool = false,
    filters: Filters = Filters { .. },
}

fn prepare(counter: MateCounter&, store: const pileup::Store&) {
    // A new region can have the same number of records but different flags.
    counter.filters_ready = false;
    counter.groups.clear();
    counter.heads.clear();
    counter.heads.resize(store.records.len, -1);
    counter.group_ends.clear();
    counter.group_ends.resize(store.records.len, -1);
    counter.overlaps.clear();
    counter.overlaps.resize(store.records.len, false);
    var capacity: i64 = 1;
    while capacity < store.records.len * 2 {
        capacity *= 2;
    }
    counter.table.clear();
    counter.table.resize(capacity, -1);
    for &record, i in store.records {
        let name = store.names[record.name..record.name + record.name_len];
        var hash: u64 = 14695981039346656037;
        for c in name {
            hash = (hash ^ c) * 1099511628211;
        }
        var slot = (hash & (capacity - 1) as u64) as i64;
        while counter.table[slot] >= 0 {
            let other = store.records[counter.table[slot]];
            if name == store.names[other.name..other.name + other.name_len] {
                break;
            }
            slot = (slot + 1) & (capacity - 1);
        }
        if counter.table[slot] < 0 {
            counter.table[slot] = i;
        }
        let group = counter.table[slot];
        counter.groups.push(ReadInfo { group: group, .. });
        // Ends are inclusive here, matching Rust's stored span. Zero-span
        // reads conservatively occupy their start for this overlap proof.
        let end = record.pos + max(record.reference_len, 1) - 1;
        if counter.group_ends[group] >= 0 {
            if record.pos <= counter.group_ends[group] {
                counter.overlaps[group] = true;
            }
            counter.group_ends[group] = max(counter.group_ends[group], end);
        } else {
            counter.group_ends[group] = end;
        }
    }
    // Mark disjoint groups once, rather than following a second array lookup
    // for every entry at every reference position.
    for &info in counter.groups {
        if !counter.overlaps[info.group] {
            info.group = -1;
        }
    }
}

// Count one reference position after resolving alignments with the same read
// name. Call counter.prepare(store) after each region load, then use this
// overload for the entries from store.columns(). All entries must come from
// that same store. mode selects the resolution rule and must be a positive
// value returned by strategy(), for example strategy("Original").
//
// ```goose
// let store = pileup::Store { .. };
// let counter = perbase::MateCounter { .. };
// counter.prepare(store);
// let row = counter.count(store, [], perbase::Filters { .. }, perbase::strategy("Original"));
// assert(row.values[0] == 0 && row.values[16] == 0);
// ```
fn count(
    counter: MateCounter&,
    store:   const pileup::Store&,
    entries: const pileup::Entry[:],
    filters: Filters,
    mode:    i64,
) -> Counts {
    // Flags and mapping quality do not change between columns. Cache their
    // filter result, but allow API callers to change filters within a region.
    // Base quality still belongs to each entry and is checked by add().
    if !counter.filters_ready ||
            counter.filters.include_flags != filters.include_flags ||
            counter.filters.exclude_flags != filters.exclude_flags ||
            counter.filters.min_mapq != filters.min_mapq {
        for record, i in store.records {
            counter.groups[i].passed =
                record.flags & filters.include_flags == filters.include_flags &&
                record.flags & filters.exclude_flags == 0 &&
                record.mapq >= filters.min_mapq;
        }
        counter.filters = filters;
        counter.filters_ready = true;
    }
    counter.links.clear();
    counter.links.resize(entries.len, -1);
    counter.touched.clear();
    var row = Counts { .. };
    for entry, i in entries {
        let info = counter.groups[entry.record_index];
        let group = info.group;
        if group < 0 {
            add(row, entry, filters, info.passed, 0, false);
            continue;
        }
        if counter.heads[group] < 0 {
            counter.touched.push(group);
        }
        counter.links[i] = counter.heads[group];
        counter.heads[group] = i;
    }
    for group in counter.touched {
        let last = counter.heads[group];
        let previous = counter.links[last];
        if previous < 0 {
            add(row, entries[last], filters, counter.groups[entries[last].record_index].passed, 0, false);
            counter.heads[group] = -1;
            continue;
        }
        if counter.links[previous] < 0 {
            // Match the existing two-element sort: compare the later entry
            // against the earlier one, then resolve the chosen ordering.
            let later_first, ignored = choose(
                store, entries[last], entries[previous],
                counter.groups[entries[last].record_index].passed,
                counter.groups[entries[previous].record_index].passed,
                mode,
            );
            let best = if later_first { last } else { previous };
            let second = if later_first { previous } else { last };
            let better, recommended = choose(
                store, entries[best], entries[second],
                counter.groups[entries[best].record_index].passed,
                counter.groups[entries[second].record_index].passed,
                mode,
            );
            add(row, entries[best], filters, counter.groups[entries[best].record_index].passed, recommended, true);
            counter.heads[group] = -1;
            continue;
        }
        counter.members.clear();
        var current = counter.heads[group];
        while current >= 0 {
            counter.members.push(current);
            current = counter.links[current];
        }
        for i in counter.members.len / 2 {
            swap(counter.members[i], counter.members[counter.members.len - i - 1]);
        }
        // Stable insertion, comparing from the back as Rust's small-slice sort
        // does. Deletion fallback makes the comparator non-transitive for some
        // groups containing more than two alignments.
        for i in 1..counter.members.len {
            let value = counter.members[i];
            var j = i;
            while j > 0 {
                let other = entries[counter.members[j - 1]];
                let better, base = choose(
                    store, entries[value], other,
                    counter.groups[entries[value].record_index].passed,
                    counter.groups[other.record_index].passed,
                    mode,
                );
                if !better {
                    break;
                }
                counter.members[j] = counter.members[j - 1];
                j -= 1;
            }
            counter.members[j] = value;
        }
        let best = counter.members[0];
        let second = if counter.members.len > 1 { counter.members[1] } else { -1 };
        var recommended: u8 = 0;
        if second >= 0 {
            let better, base = choose(
                store, entries[best], entries[second],
                counter.groups[entries[best].record_index].passed,
                counter.groups[entries[second].record_index].passed,
                mode,
            );
            recommended = base;
        }
        add(row, entries[best], filters, counter.groups[entries[best].record_index].passed, recommended, second >= 0);
        counter.heads[group] = -1;
    }
    row
}

// Decide which whole reads participate in counting, reproducing HTSlib's
// maximum-depth behavior. Fill admitted in record order and return whether
// any read was rejected. Do this before filtering or resolving mates.
//
// A read ending exactly where another starts still occupies a place for this
// decision. The first read at a new start bypasses the cap, except at reference
// ID 0, position 0. Consequently the resulting position depth may exceed limit;
// clipping each position's depth would not reproduce the same counts.
fn admit(store: const pileup::Store&, limit: i64, admitted: bool[>..]&) -> bool {
    admitted.clear();
    admitted.resize(store.records.len, true);
    if limit > 0 && store.records.len <= limit {
        return false;
    }
    var ends:     i64[>..] = [];
    var previous: i64      = -1;
    var any_rejected = false;
    for &record, i in store.records {
        let start = record.pos;
        while ends.len > 0 && ends[0] < start {
            let last = ends.pop();
            if ends.len > 0 {
                var at: i64 = 0;
                while at * 2 + 1 < ends.len {
                    var child = at * 2 + 1;
                    if child + 1 < ends.len && ends[child + 1] < ends[child] {
                        child += 1;
                    }
                    if ends[child] >= last {
                        break;
                    }
                    ends[at] = ends[child];
                    at = child;
                }
                ends[at] = last;
            }
        }
        let at_cursor = start == previous || (previous < 0 && record.tid == 0 && start == 0);
        if at_cursor && ends.len >= limit {
            admitted[i] = false;
            any_rejected = true;
            continue;
        }
        let end = start + record.reference_len;
        if end > start || !at_cursor {
            var at = ends.len;
            ends.push(end);
            while at > 0 && ends[(at - 1) / 2] > end {
                ends[at] = ends[(at - 1) / 2];
                at = (at - 1) / 2;
            }
            ends[at] = end;
        }
        previous = start;
    }
    any_rejected
}

// Append one tab-separated output row and a newline. Quote reference names
// containing tabs, newlines or quotes, doubling any quote inside the name.
// position must already use the desired zero- or one-based coordinates.
// reference is the reference-base character; zero omits the REF_BASE column.
// This only appends bytes to output; the caller writes or flushes them.
//
// ```goose
// var text: u8[>..] = [];
// perbase::append_row(text, "chr1", 7, 0, perbase::Counts { .. }, false);
// assert(text[..7] == "chr1\t7\t" && text[text.len - 6..] == "false\n");
// ```
fn append_row(
    output: u8[>..]&, name: const u8[:], position: i64, reference: u8,
    row:    Counts, near_max: bool,
) {
    var quoted = false;
    for c in name {
        if c == '\t' || c == '\r' || c == '\n' || c == '"' {
            quoted = true;
        }
    }
    if quoted {
        output.push('"');
        for c in name {
            if c == '"' {
                output.push('"');
            }
            output.push(c);
        }
        output.push('"');
    } else {
        output.append(name);
    }
    format(output, "\t", position);
    if reference != 0 {
        output.push('\t');
        output.push(reference);
    }
    for value in row.values {
        format(output, "\t", value);
    }
    output.append(if near_max { "\ttrue\n" } else { "\tfalse\n" });
}
