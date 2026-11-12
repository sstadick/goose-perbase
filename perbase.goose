// Perbase-compatible base-depth aggregation over Seqair-style raw columns.
namespace perbase;
import std;
import seqair from "github.com/sstadick/goose-seqair@v0.1.0";
import seqair.bam;
import seqair.pileup;

// Strategies use Perbase's case-insensitive CLI names. Zero disables resolving.
fn strategy(name: const u8[:]) -> i64 {
    let names = ["original", "basequalmapqualfirstinpair", "basequalmapqualiupac",
        "basequalmapqualn", "mapqualbasequalfirstinpair", "mapqualbasequaliupac",
        "mapqualbasequaln", "iupac", "n"];
    for candidate, i in names {
        if candidate.len != name.len { continue; }
        var equal = true;
        for c, j in name { if (if c >= 'A' && c <= 'Z' { c + 32 } else { c }) != candidate[j] { equal = false; } }
        if equal { return i + 1; }
    }
    -1
}

// Apply read filters at the column, not at loading: rejected reads contribute
// to FAIL. Defaults match Perbase, not goose-seqair's default exclusion mask.
// A base below min_base_quality counts as N without reducing depth.
//
// ```goose
// let options = perbase::Filters { .. };
// assert(options.exclude_flags == 0 && options.min_base_quality == 0);
// ```
struct Filters {
    include_flags: u16 = 0, exclude_flags: u16 = 0,
    min_mapq: u8 = 0, min_base_quality: u8 = 0,
}

// The 17 numeric columns after POS (or REF_BASE), in Perbase TSV order.
// Deletions contribute to DEPTH; reference skips do not. INS counts events,
// not inserted bases. IUPAC ambiguity codes in reads count as N, like Perbase.
// The final counter counts passing query-name groups with multiple alignments.
//
// ```goose
// let row = perbase::Counts { .. };
// assert(row.values[0] == 0 && row.values[16] == 0);
// ```
struct Counts { values: i64[17] = [0; 17] }

// Count a raw column using its store's metadata. Entries borrow the store;
// the returned counters own no input views. Missing quality 255 is compared
// numerically, matching Perbase's behavior rather than inventing a cutoff.
//
// ```goose
// let store = pileup::Store { .. };
// let row = perbase::count(store, [], perbase::Filters { .. });
// assert(row.values[0] == 0);
// ```
fn count(store: const pileup::Store&, entries: const pileup::Entry[:], filters: Filters) -> Counts {
    var row = Counts { .. };
    for &entry in entries {
        add(row, store, entry, filters, 0, false);
    }
    row
}

fn passes(store: const pileup::Store&, entry: pileup::Entry, filters: Filters) -> bool {
    let f = store.records[entry.record_index];
    (f.flags & filters.include_flags) == filters.include_flags && f.flags & filters.exclude_flags == 0 &&
        f.mapq >= filters.min_mapq
}

fn add(row: Counts&, store: const pileup::Store&, entry: pileup::Entry, filters: Filters, recommended: u8, resolved: bool) {
    if !passes(store, entry, filters) { row.values[15] += 1; return; }
    if resolved { row.values[16] += 1; }
    if entry.op == pileup::Op.RefSkip { row.values[14] += 1; return; }
    row.values[0] += 1;
    if entry.op == pileup::Op.Deletion { row.values[13] += 1; return; }
    var base = if recommended != 0 { recommended } else { entry.base };
    if entry.quality < filters.min_base_quality { base = 'N'; }
    var column: i64 = 5;
    if base == 'A' { column = 1; } else if base == 'C' { column = 2; }
    else if base == 'G' { column = 3; } else if base == 'T' || base == 'U' { column = 4; }
    else if recommended != 0 {
        for c, i in "RYSWKM" { if c == base { column = 6 + i; } }
    }
    row.values[column] += 1;
    if entry.indel > 0 { row.values[12] += 1; }
}

fn canonical(base: u8) -> u8 {
    if base == 'A' || base == 'C' || base == 'G' || base == 'T' { base }
    else if base == 'U' { 'T' } else { 'N' }
}

fn choose(store: const pileup::Store&, a: pileup::Entry, b: pileup::Entry, filters: Filters, mode: i64) -> bool, u8 {
    let ap = passes(store, a, filters); let bp = passes(store, b, filters);
    if ap != bp { return ap, 0; }
    if !ap { return true, 0; }
    let af = store.records[a.record_index]; let bf = store.records[b.record_index];
    let fallback = mode == 1 || a.qpos < 0 || b.qpos < 0;
    if fallback {
        if af.mapq != bf.mapq { return af.mapq > bf.mapq, 0; }
    } else if mode < 8 {
        let aq = if mode <= 4 { a.quality } else { af.mapq };
        let bq = if mode <= 4 { b.quality } else { bf.mapq };
        if aq != bq { return aq > bq, 0; }
        let ar = if mode <= 4 { af.mapq } else { a.quality };
        let br = if mode <= 4 { bf.mapq } else { b.quality };
        if ar != br { return ar > br, 0; }
    }
    if fallback || mode == 2 || mode == 5 { return af.flags & 64 != 0 || bf.flags & 64 == 0, 0; }
    let ab = canonical(a.base); let bb = canonical(b.base);
    if ab == bb { return true, ab; }
    if mode == 4 || mode == 7 || mode == 9 { return true, 'N'; }
    let lo = min(ab, bb); let hi = max(ab, bb);
    var base: u8 = 'N';
    if lo == 'A' && hi == 'G' { base = 'R'; }
    else if lo == 'C' && hi == 'T' { base = 'Y'; }
    else if lo == 'C' && hi == 'G' { base = 'S'; }
    else if lo == 'A' && hi == 'T' { base = 'W'; }
    else if lo == 'G' && hi == 'T' { base = 'K'; }
    else if lo == 'A' && hi == 'C' { base = 'M'; }
    return true, base;
}

// Query-name grouping is built once per region, including secondary/unpaired
// records: Perbase groups by name, not solely by Seqair's eligible mate links.
root struct MateCounter {
    groups: i64[>..] = [], heads: i64[>..] = [], links: i64[>..] = [], touched: i64[>..] = [],
    members: i64[>..] = [],
}

fn prepare(counter: MateCounter&, store: const pileup::Store&) {
    counter.groups.clear(); counter.heads.clear();
    counter.heads.resize(store.records.len, -1);
    var capacity: i64 = 1;
    while capacity < store.records.len * 2 { capacity *= 2; }
    var table: i64[>..] = []; table.resize(capacity, -1);
    for &record, i in store.records {
        let name = store.names[record.name..record.name + record.name_len];
        var hash: u64 = 14695981039346656037;
        for c in name { hash = (hash ^ c) * 1099511628211; }
        var slot = (hash & (capacity - 1) as u64) as i64;
        while table[slot] >= 0 {
            let other = store.records[table[slot]];
            if name == store.names[other.name..other.name + other.name_len] { break; }
            slot = (slot + 1) & (capacity - 1);
        }
        if table[slot] < 0 { table[slot] = i; }
        counter.groups.push(table[slot]);
    }
}

// Resolve query-name groups before counting. Call counter.prepare(store) once
// after each region load, then reuse it across that store's columns. Mode is a
// positive value returned by strategy(); input entries must belong to store.
//
// ```goose
// let store = pileup::Store { .. };
// let counter = perbase::MateCounter { .. }; counter.prepare(store);
// let row = counter.count(store, [], perbase::Filters { .. }, perbase::strategy("Original"));
// assert(row.values[0] == 0 && row.values[16] == 0);
// ```
fn count(counter: MateCounter&, store: const pileup::Store&, entries: const pileup::Entry[:], filters: Filters, mode: i64) -> Counts {
    counter.links.clear(); counter.links.resize(entries.len, -1); counter.touched.clear();
    for entry, i in entries {
        let group = counter.groups[entry.record_index];
        if counter.heads[group] < 0 { counter.touched.push(group); }
        counter.links[i] = counter.heads[group]; counter.heads[group] = i;
    }
    var row = Counts { .. };
    for group in counter.touched {
        counter.members.clear(); var current = counter.heads[group];
        while current >= 0 { counter.members.push(current); current = counter.links[current]; }
        for i in counter.members.len / 2 { swap(counter.members[i], counter.members[counter.members.len - i - 1]); }
        // Stable insertion, comparing from the back as Rust's small-slice sort
        // does. Deletion fallback makes the comparator non-transitive for some
        // groups containing more than two alignments.
        for i in 1..counter.members.len {
            let value = counter.members[i]; var j = i;
            while j > 0 {
                let better, base = choose(store, entries[value], entries[counter.members[j - 1]], filters, mode);
                if !better { break; }
                counter.members[j] = counter.members[j - 1]; j -= 1;
            }
            counter.members[j] = value;
        }
        let best = counter.members[0]; let second = if counter.members.len > 1 { counter.members[1] } else { -1 };
        var recommended: u8 = 0;
        if second >= 0 { let better, base = choose(store, entries[best], entries[second], filters, mode); recommended = base; }
        add(row, store, entries[best], filters, recommended, second >= 0);
        counter.heads[group] = -1;
    }
    row
}

// Reproduce HTSlib's read-admission cap, not a per-column truncation. End==start
// entries remain counted during admission at that start. A new start's first
// record bypasses the cap, except the initial cursor (TID 0, position 0).
fn admit(store: const pileup::Store&, limit: i64, admitted: bool[>..]&) -> bool {
    admitted.clear(); admitted.resize(store.records.len, true);
    if limit > 0 && store.records.len <= limit { return false; }
    var ends: i64[>..] = []; var previous: i64 = -1; var any_rejected = false;
    for &record, i in store.records {
        let start = record.pos;
        while ends.len > 0 && ends[0] < start {
            let last = ends.pop();
            if ends.len > 0 {
                var at: i64 = 0;
                while at * 2 + 1 < ends.len {
                    var child = at * 2 + 1;
                    if child + 1 < ends.len && ends[child + 1] < ends[child] { child += 1; }
                    if ends[child] >= last { break; }
                    ends[at] = ends[child]; at = child;
                }
                ends[at] = last;
            }
        }
        let at_cursor = start == previous || (previous < 0 && record.tid == 0 && start == 0);
        if at_cursor && ends.len >= limit { admitted[i] = false; any_rejected = true; continue; }
        let end = start + record.reference_len;
        if end > start || !at_cursor {
            var at = ends.len; ends.push(end);
            while at > 0 && ends[(at - 1) / 2] > end { ends[at] = ends[(at - 1) / 2]; at = (at - 1) / 2; }
            ends[at] = end;
        }
        previous = start;
    }
    any_rejected
}

// Append a TSV row, quoting reference names as CSV with a tab delimiter.
// Position is already in the requested coordinate base; reference=0 omits
// REF_BASE. The caller controls buffering and I/O.
//
// ```goose
// var text: u8[>..] = [];
// perbase::append_row(text, "chr1", 7, 0, perbase::Counts { .. }, false);
// assert(text[..7] == "chr1\t7\t" && text[text.len - 6..] == "false\n");
// ```
fn append_row(output: u8[>..]&, name: const u8[:], position: i64, reference: u8,
              row: Counts, near_max: bool) {
    var quoted = false;
    for c in name { if c == '\t' || c == '\r' || c == '\n' || c == '"' { quoted = true; } }
    if quoted {
        output.push('"');
        for c in name { if c == '"' { output.push('"'); } output.push(c); }
        output.push('"');
    } else { output.append(name); }
    format(output, "\t", position);
    if reference != 0 { output.push('\t'); output.push(reference); }
    for value in row.values { format(output, "\t", value); }
    output.append(if near_max { "\ttrue\n" } else { "\tfalse\n" });
}
