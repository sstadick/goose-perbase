// Region selection needs only BED outer intervals or VCF/BCF record positions.
// BCF genotypes/INFO are skipped, not decoded into an application object model.
namespace regions;
import std;
import seqair from "github.com/sstadick/goose-seqair@v0.1.0";
import seqair.bam;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

struct Interval { tid: i32, begin: i64, end: i64, ordinal: i64 = 0 }
fn check(ok: bool) { if !ok { write_stderr("goose-perbase: invalid interval input or unknown reference\n"); exit(1); } }

root struct Input {
    file: io::FileHandle = io::FileHandle { fd: -1 }, zipped: bgzf::Reader,
    compressed: bool = false, buffer: u8[65536], begin: i64 = 0, end: i64 = 0,
}

fn open(input: Input&, path: const u8[:]) {
    let file, ok = io::open_read(path); check(ok); input.file.fd = file.fd; file.fd = -1;
    let n = io::read(input.file, input.buffer[..2]); check(n >= 0);
    input.compressed = n == 2 && input.buffer[0] == 31 && input.buffer[1] == 139;
    if input.compressed { io::close(input.file); check(input.zipped.open(path) == bgzf::Status.Ready); }
    else { check(bgzf::seek_raw(input.file.fd, 0) == 0); }
}

fn close(input: Input&) { io::close(input.file); input.zipped.close(); }
fn fill(input: Input&) -> bool {
    if input.begin < input.end { return true; }
    input.begin = 0;
    if input.compressed {
        let n, status = input.zipped.read_some(input.buffer);
        check(status == bgzf::Status.Ready || status == bgzf::Status.End); input.end = n;
    } else { input.end = io::read(input.file, input.buffer); check(input.end >= 0); }
    input.end > 0
}
fn exact(input: Input&, bytes: u8[:]) {
    var at: i64 = 0;
    while at < bytes.len {
        check(input.fill()); let n = min(bytes.len - at, input.end - input.begin);
        bytes[at..at + n].move_from(input.buffer[input.begin..input.begin + n]);
        at += n; input.begin += n;
    }
}
fn skip(input: Input&, count: i64) {
    var left = count;
    while left > 0 { check(input.fill()); let n = min(left, input.end - input.begin); input.begin += n; left -= n; }
}
fn line(input: Input&, bytes: u8[>..]&) -> bool {
    bytes.clear(); var consumed = false;
    while input.fill() {
        consumed = true;
        let start = input.begin;
        while input.begin < input.end && input.buffer[input.begin] != '\n' { input.begin += 1; }
        bytes.append(input.buffer[start..input.begin]);
        if input.begin < input.end { input.begin += 1; break; }
    }
    if bytes.len > 0 && bytes[bytes.len - 1] == '\r' { bytes.resize(bytes.len - 1); }
    consumed
}
fn number(bytes: const u8[:]) -> i64 {
    check(bytes.len > 0); var n: i64 = 0;
    for c in bytes { check(c >= '0' && c <= '9' && n <= (4294967295 - (c as i64 - '0')) / 10); n = n * 10 + c - '0'; }
    n
}
fn add<R>(intervals: Interval[>..]&, reader: const R&, name: const u8[:], begin: i64, end: i64) {
    let tid = reader.reference_id(name); check(tid >= 0 && end >= begin);
    if begin < end && begin < reader.references[tid as i64].length {
        intervals.push(Interval { tid: tid, begin: begin, end: min(end, reader.references[tid as i64].length), ordinal: intervals.len });
    }
}

fn read<R>(path: const u8[:], variants: bool, reader: const R&, intervals: Interval[>..]&) {
    let input = Input { .. }; input.open(path); defer input.close();
    if !input.fill() { return; }
    if variants && input.buffer[input.begin] == 'B' {
        var prefix: u8[9] = [0; 9]; input.exact(prefix);
        check(prefix[..5] == "BCF\x02\x02");
        let len = bgzf::le32(prefix, 5) as i64; check(len > 0 && len <= 16777216);
        var text: u8[>..] = []; text.resize(len, 0); input.exact(text);
        var tids: i32[>..] = [];
        bstr::split_str(text, "\n") { row =>
            if starts_with(row, "##contig=<") {
                var tid: i32 = -1; var idx = tids.len;
                let stop = if row.len > 0 && row[row.len - 1] == '\r' { row.len - 2 } else { row.len - 1 };
                bstr::split_str(row[10..stop], ",") { field =>
                    if starts_with(field, "ID=") { tid = reader.reference_id(field[3..]); }
                    if starts_with(field, "IDX=") { idx = number(field[4..]); }
                };
                check(idx <= 1048576);
                if tids.len <= idx { tids.resize(idx + 1, -1); }
                tids[idx] = tid;
            }
        };
        var lengths: u8[8] = [0; 8]; var fixed: u8[8] = [0; 8];
        while input.fill() {
            input.exact(lengths); let shared = bgzf::le32(lengths, 0) as i64; let individual = bgzf::le32(lengths, 4) as i64;
            check(shared >= 24); input.exact(fixed);
            let rid = bgzf::le32(fixed, 0) as i64; let pos = bgzf::le32(fixed, 4) as i64;
            check(rid < tids.len && tids[rid] >= 0 && pos < 2147483647);
            let ref = reader.references[tids[rid] as i64];
            add(intervals, reader, reader.names[ref.name..ref.name + ref.name_len], pos, pos + 1);
            input.skip(shared - 8 + individual);
        }
        return;
    }
    var text: u8[>..] = [];
    while input.line(text) {
        if text.len == 0 || text[0] == '#' || starts_with(text, "track ") || starts_with(text, "browser ") { continue; }
        var fields: const u8[:][>..] = [];
        bstr::split_str(text, "\t") { field => fields.push(field); };
        check(fields.len >= (if variants { 8 } else { 3 }));
        let start = number(fields[1]);
        if variants { check(start > 0); add(intervals, reader, fields[0], start - 1, start); }
        else { add(intervals, reader, fields[0], start, number(fields[2])); }
    }
}

fn merge(intervals: Interval[>..]&, merge_overlaps: bool) {
    sort(intervals) { a, b => a.tid < b.tid || (a.tid == b.tid && (a.begin < b.begin || (a.begin == b.begin && (a.end < b.end || (a.end == b.end && a.ordinal < b.ordinal))))) };
    if !merge_overlaps { return; }
    var kept: i64 = 0;
    for interval in intervals {
        if kept > 0 && intervals[kept - 1].tid == interval.tid && intervals[kept - 1].end >= interval.begin {
            intervals[kept - 1].end = max(intervals[kept - 1].end, interval.end);
        } else { intervals[kept] = interval; kept += 1; }
    }
    intervals.resize(kept);
}
