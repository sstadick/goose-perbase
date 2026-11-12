// Indexed input over worker-owned handles; counting consumes borrowed BAM records.
namespace alignment;
import std;
import os;
import seqair from "github.com/sstadick/goose-seqair@v0.2.1";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import seqair.cram;
import seqair.cram_index;
import seqair.sam;
import bio from "github.com/sstadick/goose-bio@v0.4.1";
import bio.faidx;
import io from "github.com/sstadick/goose-io@v0.3.0";

struct Reference { reader: faidx::Reader& }
fn fetch(source: Reference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool {
    source.reader.fetch(name, begin, end, output) == faidx::Status.Ok
}

root struct Header {
    is_cram: bool = false,
    references: bam::Reference[>..] = [], names: u8[>..] = [], index_path: u8[>..] = [],
    status: bam::Status = bam::Status.Closed,
}

root struct BamCursor {
    chunks: bam_index::Chunk[>..] = [], chunk: i64 = 0, active: bool = false,
    tid: i32 = -1, begin: i64 = 0, end: i64 = 0, last_tid: i32 = -1, last_pos: i64 = -1,
    status: bam::Status = bam::Status.Closed,
}

// Pass owners directly when opening/closing: the pinned Goose compiler merges
// shrink effects through stored root references with borrowed CLI paths.
fn open(reader: Header&, bam_reader: bam::Reader&, index: bam_index::Index&,
        cram_reader: cram::Reader&, crai: cram_index::Index&, reference: faidx::Reader&,
        path: const u8[:], index_path: const u8[:] = "", reference_path: const u8[:] = "") -> bam::Status {
    bam_reader.close(); cram_reader.close(); reference.close();
    reader.references.clear(); reader.names.clear(); reader.index_path.clear();
    let file, ok = io::open_read(path); if !ok { return bam::Status.IoError; }
    var magic: u8[4] = [0; 4]; let n = io::read(file, magic); io::close(file);
    reader.is_cram = n == 4 && magic == "CRAM";
    if reader.is_cram {
        reader.status = cram_reader.open(path); if reader.status != bam::Status.Record { return reader.status; }
        reader.references.append(cram_reader.references); reader.names.append(cram_reader.names);
    } else {
        reader.status = bam_reader.open(path); if reader.status != bam::Status.Record { return reader.status; }
        reader.references.append(bam_reader.references); reader.names.append(bam_reader.names);
    }
    if index_path.len > 0 { reader.index_path.append(index_path); }
    else if reader.is_cram {
        format(reader.index_path, path, ".crai");
        if !file_exists(reader.index_path) && ends_with(path, ".cram") {
            reader.index_path.clear(); format(reader.index_path, path[..path.len - 5], ".crai");
        }
    }
    else {
        format(reader.index_path, path, ".bai");
        if !file_exists(reader.index_path) {
            reader.index_path.clear(); format(reader.index_path, path, ".csi");
            if !file_exists(reader.index_path) && path.len > 4 {
                reader.index_path.clear(); format(reader.index_path, path[..path.len - 4], ".bai");
            }
        }
    }
    if reader.is_cram {
        reader.status = crai.open(reader.index_path);
        for entry in crai.entries { if entry.tid >= reader.references.len { reader.status = bam::Status.Invalid; } }
    } else {
        reader.status = index.open(reader.index_path);
        if reader.status == bam::Status.Record && index.references.len != reader.references.len { reader.status = bam::Status.Invalid; }
    }
    if reader.status == bam::Status.Record && reference_path.len > 0 && reference.open(reference_path) != faidx::Status.Ok { reader.status = bam::Status.Invalid; }
    reader.status
}

fn reference_id(reader: const Header&, name: const u8[:]) -> i32 { sam::find_reference(reader.references, reader.names, name) }

fn fetch(reader: bam::Reader&, index: bam_index::Index&, cursor: BamCursor&, tid: i32, begin: i64, end: i64) -> bam::Status {
    cursor.tid = tid; cursor.begin = begin; cursor.end = end; cursor.chunk = 0; cursor.active = false;
    cursor.last_tid = -1; cursor.last_pos = -1;
    cursor.status = index.query(tid, begin, end, cursor.chunks); cursor.status
}

fn read(reader: bam::Reader&, cursor: BamCursor&, source: Reference&) -> bam::Status {
    if cursor.status != bam::Status.Record { return cursor.status; }
    while cursor.chunk < cursor.chunks.len {
        let chunk = cursor.chunks[cursor.chunk];
        if !cursor.active {
            cursor.status = reader.seek(chunk.begin);
            if cursor.status != bam::Status.Record { return cursor.status; } cursor.active = true;
        }
        if reader.tell() >= chunk.end { cursor.chunk += 1; cursor.active = false; continue; }
        cursor.status = reader.read();
        if cursor.status == bam::Status.End { cursor.status = bam::Status.Truncated; }
        if cursor.status != bam::Status.Record { return cursor.status; }
        let f = reader.fields;
        if f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < cursor.last_tid || (f.tid == cursor.last_tid && f.pos < cursor.last_pos) { cursor.status = bam::Status.Unsorted; return cursor.status; }
        cursor.last_tid = f.tid; cursor.last_pos = f.pos;
        if f.tid > cursor.tid || (f.tid == cursor.tid && f.pos >= cursor.end) { cursor.chunk += 1; cursor.active = false; continue; }
        if f.tid == cursor.tid && f.end > cursor.begin { return cursor.status; }
    }
    cursor.status = bam::Status.End; cursor.status
}

fn fetch(reader: cram::Reader&, index: cram_index::Index&, cursor: cram_index::Cursor&, tid: i32, begin: i64, end: i64) -> bam::Status {
    cursor.fetch(index, reader, tid, begin, end)
}

fn read(reader: cram::Reader&, cursor: cram_index::Cursor&, source: Reference&) -> bam::Status {
    cursor.read(reader, source)
}

fn load(store: pileup::Store&, reader: cram::Reader&, index: cram_index::Index&, cursor: cram_index::Cursor&,
        source: Reference&, tid: i32, begin: i64, end: i64) -> bam::Status {
    store.clear(); store.tid = tid; store.begin = begin; store.end = end;
    var status = cursor.fetch(index, reader, tid, begin, end);
    while status == bam::Status.Record {
        status = cursor.read(reader, source);
        if status != bam::Status.Record { break; }
        let record = reader.record();
        if record.fields.flags & 4 == 0 { status = pileup::retain(store, record); }
    }
    if status != bam::Status.End { store.clear(); return status; }
    store.link_mates(); bam::Status.Record
}
