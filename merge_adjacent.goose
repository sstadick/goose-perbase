// Streaming Perbase TSV merging; memory is bounded by a record and output batch.
namespace merge_adjacent;
import std;
import regions;
import output;
import opts from "github.com/sstadick/goose-opts@v0.1.0";
import io from "github.com/sstadick/goose-io@v0.3.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

struct Options {
    in_file: opts::Field<const u8[:]> = opts::Field<const u8[:]> { positional: true, help: "Sorted TSV; default stdin (-)" },
    output: opts::Field<const u8[:]> = opts::Field<const u8[:]> { short: "o", help: "Output file; default stdout" },
    no_header: opts::Field<bool> = opts::Field<bool> { long: "no-header", short: "n", help: "Input columns are REF, POS, END, DEPTH without a header" },
    bgzip: opts::Field<bool> = opts::Field<bool> { short: "Z", help: "BGZF-compress output" },
    compression_threads: opts::Field<i64> = opts::Field<i64> { value: hardware_threads(), long: "compression-threads", short: "T", help: "Output compression workers" },
    compression_level: opts::Field<i64> = opts::Field<i64> { value: 2, long: "compression-level", short: "L", help: "Compression level 0..9" },
}

fn check(ok: bool, message: const u8[:]) {
    if !ok { write_stderr(str("goose-perbase merge-adjacent: ", message, "\n")); exit(1); }
}

root struct Table { bytes: u8[>..] = [], ends: i64[>..] = [], width: i64 = -1, skip_lf: bool = false }

// CSV quoting with tab delimiters, including escaped quotes and multiline names.
// Like csv::Reader, ignore blank lines and require a consistent field count.
fn next(table: Table&, input: regions::Input&) -> bool {
    table.bytes.clear(); table.ends.clear();
    var quoted = false; var quote_end = false; var field_start = true; var started = false;
    loop {
        if !input.fill() { if !started { return false; } break; }
        let c = input.buffer[input.begin]; input.begin += 1;
        if table.skip_lf { table.skip_lf = false; if c == '\n' { continue; } }
        if quoted {
            if c == '"' { quoted = false; quote_end = true; }
            else { table.bytes.push(c); }
            continue;
        }
        if quote_end && c == '"' { table.bytes.push(c); quoted = true; quote_end = false; continue; }
        quote_end = false;
        if c == '\r' || c == '\n' {
            table.skip_lf = c == '\r';
            if !started { continue; }
            break;
        }
        started = true;
        if c == '\t' { table.ends.push(table.bytes.len); field_start = true; }
        else {
            if c == '"' && field_start { quoted = true; }
            else { table.bytes.push(c); }
            field_start = false;
        }
    }
    table.ends.push(table.bytes.len);
    if table.width < 0 { table.width = table.ends.len; }
    check(table.ends.len == table.width, "inconsistent TSV field count");
    true
}

fn field(table: const Table&, column: i64) -> const u8[:] {
    if column < 0 { return ""; }
    table.bytes[(if column == 0 { 0 } else { table.ends[column - 1] })..table.ends[column]]
}

fn number(bytes: const u8[:]) -> u64 {
    let start = if bytes.len > 0 && bytes[0] == '+' { 1 } else { 0 };
    check(bytes.len > start, "expected an unsigned integer"); var n: u64 = 0;
    for c in bytes[start..] {
        check(c >= '0' && c <= '9', "expected an unsigned integer");
        let digit = (c - '0') as u64;
        check(n <= (18446744073709551615 as u64 - digit) / 10, "integer overflow"); n = n * 10 + digit;
    }
    n
}

fn emit(text: u8[>..]&, name: const u8[:], begin: u64, end: u64, depth: u64) {
    var quoted = false;
    for c in name { if c == '\t' || c == '\r' || c == '\n' || c == '"' { quoted = true; } }
    if quoted { text.push('"'); for c in name { if c == '"' { text.push('"'); } text.push(c); } text.push('"'); }
    else { text.append(name); }
    format(text, "\t", begin, "\t", end, "\t", depth, "\n");
}

fn run(options: const Options&) {
    check(options.compression_threads.value > 0, "compression threads must be positive");
    check(options.compression_level.value >= 0 && options.compression_level.value <= 9, "invalid compression level");
    let path = options.in_file.value;
    check(options.output.value.len == 0 || options.output.value == "-" || options.output.value != path, "output would overwrite input");
    let input = regions::Input { .. }; defer input.close();
    if path.len == 0 || path == "-" { input.file.fd = 0; }
    else {
        input.compressed = ends_with(path, ".gz") || ends_with(path, ".gzip") || ends_with(path, ".bgzf");
        if input.compressed { check(input.zipped.open(path) == bgzf::Status.Ready, "cannot open BGZF input"); }
        else { let file, ok = io::open_read(path); check(ok, "cannot open input"); input.file.fd = file.fd; file.fd = -1; }
    }
    // csv::Reader consumes a UTF-8 BOM only at the beginning of the stream.
    input.fill();
    while input.end > 0 && input.end < 3 {
        var n: i64 = 0;
        if input.compressed {
            let got, status = input.zipped.read_some(input.buffer[input.end..3]);
            check(status == bgzf::Status.Ready || status == bgzf::Status.End, "cannot read BGZF input"); n = got;
        } else { n = io::read(input.file, input.buffer[input.end..3]); check(n >= 0, "cannot read input"); }
        if n == 0 { break; } input.end += n;
    }
    if input.end >= 3 && input.buffer[..3] == "\xef\xbb\xbf" { input.begin = 3; }
    let sink = output::Sink { .. };
    sink.open(options.output.value, options.bgzip.value, min(options.compression_threads.value, hardware_threads()), options.compression_level.value as i32);
    let table = Table { .. };
    var columns: i64[4] = [0, 1, 2, 3];
    if !options.no_header.value && table.next(input) {
        columns = [-1, -1, -1, -1];
        for i in table.ends.len {
            let name = table.field(i);
            let slot = if name == "REF" || name == "chrom" { 0 } else if name == "POS" || name == "chromStart" { 1 }
                else if name == "END" || name == "chromEnd" { 2 } else if name == "DEPTH" || name == "COV" { 3 } else { -1 };
            if slot >= 0 { check(columns[slot] < 0, "duplicate header field"); columns[slot] = i; }
        }
    }
    var name: u8[>..] = []; var text: u8[>..] = [];
    var begin: u64 = 0; var end: u64 = 0; var depth: u64 = 0; var active = false;
    while table.next(input) {
        check(bstr::is_utf8(table.bytes), "invalid UTF-8 input");
        for i in [0, 1, 3] { check(columns[i] >= 0 && columns[i] < table.ends.len, "missing REF, POS or DEPTH field"); }
        check(columns[2] < table.ends.len, "missing END field");
        let ref = table.field(columns[0]); let start = number(table.field(columns[1]));
        let stop_text = table.field(columns[2]); let stop = if stop_text.len > 0 { number(stop_text) } else { start + 1 };
        let count = number(table.field(columns[3]));
        if active && name == ref && depth == count && end >= start {
            // Upstream replaces the end, even for a nested shorter interval.
            end = stop;
        } else {
            if active { emit(text, name, begin, end, depth); }
            else { text.append("REF\tPOS\tEND\tDEPTH\n"); }
            name.clear(); name.append(ref); begin = start; end = stop; depth = count; active = true;
        }
        if text.len >= 65536 { sink.write(text); text.clear(); }
    }
    if active { emit(text, name, begin, end, depth); }
    sink.write(text); sink.close();
}
