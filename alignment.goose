// Indexed input over worker-owned handles; counting consumes borrowed BAM records.
namespace alignment;
import std;
import os;
import seqair from "github.com/sstadick/goose-seqair@v0.1.0";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import bio from "github.com/sstadick/goose-bio@v0.3.1";
import bio.faidx;

root struct Reader {
    bam_reader: bam::Reader&, index: bam_index::Index&,
    reference: faidx::Reader&,
    references: bam::Reference[>..] = [], names: u8[>..] = [], index_path: u8[>..] = [],
    chunks: bam_index::Chunk[>..] = [], chunk: i64 = 0, active: bool = false,
    tid: i32 = -1, begin: i64 = 0, end: i64 = 0, last_tid: i32 = -1, last_pos: i64 = -1,
    status: bam::Status = bam::Status.Closed,
}

fn open(reader: Reader&, bam_reader: bam::Reader&, index: bam_index::Index&, reference: faidx::Reader&,
        path: const u8[:], index_path: const u8[:] = "", reference_path: const u8[:] = "") -> bam::Status {
    bam_reader.close(); reference.close(); reader.references.clear(); reader.names.clear(); reader.index_path.clear();
    reader.status = bam_reader.open(path); if reader.status != bam::Status.Record { return reader.status; }
    reader.references.append(bam_reader.references); reader.names.append(bam_reader.names);
    if index_path.len > 0 { reader.index_path.append(index_path); }
    else {
        format(reader.index_path, path, ".bai");
        if !file_exists(reader.index_path) {
            reader.index_path.clear(); format(reader.index_path, path, ".csi");
            if !file_exists(reader.index_path) && path.len > 4 {
                reader.index_path.clear(); format(reader.index_path, path[..path.len - 4], ".bai");
            }
        }
    }
    reader.status = index.open(reader.index_path);
    if reader.status == bam::Status.Record && index.references.len != reader.references.len { reader.status = bam::Status.Invalid; }
    if reader.status == bam::Status.Record && reference_path.len > 0 && reference.open(reference_path) != faidx::Status.Ok { reader.status = bam::Status.Invalid; }
    reader.status
}

fn reference_id(reader: const Reader&, name: const u8[:]) -> i32 { reader.bam_reader.reference_id(name) }

fn fetch(reader: Reader&, tid: i32, begin: i64, end: i64) -> bam::Status {
    reader.tid = tid; reader.begin = begin; reader.end = end; reader.chunk = 0; reader.active = false;
    reader.last_tid = -1; reader.last_pos = -1;
    reader.status = reader.index.query(tid, begin, end, reader.chunks); reader.status
}

fn read(reader: Reader&) -> bam::Status {
    if reader.status != bam::Status.Record { return reader.status; }
    while reader.chunk < reader.chunks.len {
        let chunk = reader.chunks[reader.chunk];
        if !reader.active {
            reader.status = reader.bam_reader.seek(chunk.begin);
            if reader.status != bam::Status.Record { return reader.status; } reader.active = true;
        }
        if reader.bam_reader.tell() >= chunk.end { reader.chunk += 1; reader.active = false; continue; }
        reader.status = reader.bam_reader.read();
        if reader.status == bam::Status.End { reader.status = bam::Status.Truncated; }
        if reader.status != bam::Status.Record { return reader.status; }
        let f = reader.bam_reader.fields;
        if f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < reader.last_tid || (f.tid == reader.last_tid && f.pos < reader.last_pos) { reader.status = bam::Status.Unsorted; return reader.status; }
        reader.last_tid = f.tid; reader.last_pos = f.pos;
        if f.tid > reader.tid || (f.tid == reader.tid && f.pos >= reader.end) { reader.chunk += 1; reader.active = false; continue; }
        if f.tid == reader.tid && f.end > reader.begin { return reader.status; }
    }
    reader.status = bam::Status.End; reader.status
}

fn record(reader: const Reader&) -> bam::Record { reader.bam_reader.record() }

fn load(reader: Reader&, store: pileup::Store&, tid: i32, begin: i64, end: i64) -> bam::Status {
    store.load_indexed(reader.bam_reader, reader.index, tid, begin, end,
        pileup::Filter { exclude_flags: 0, min_mapq: 0, retain_zero_span: true })
}
