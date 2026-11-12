// Ordered, bounded BGZF compression using persistent Goose workers.
namespace output;
import std;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";

struct Block { id: i64, data: u8[] }
struct Encoded { id: i64, data: u8[] }

fn check(ok: bool) { if !ok { write_stderr("goose-perbase: output/compression failure\n"); exit(1); } }

thread_fn compress(level: i32) {
    let encoder, ok = deflate::compressor(level); check(ok); defer encoder.close();
    var bytes: u8[65536] = [0; 65536];
    loop {
        let message = qget<Block>();
        if message.id < 0 { break; }
        let n = encoder.compress(deflate::Format.Raw, message.data, bytes[18..65528]);
        check(n > 0);
        bytes[..18].move_from("\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0\0\0");
        bgzf::put16(bytes, 16, (n + 25) as u32);
        bgzf::put32(bytes, n + 18, deflate::crc32(message.data));
        bgzf::put32(bytes, n + 22, message.data.len as u32);
        check(qput(Encoded { id: message.id, data: bytes[..n + 26] }));
    }
}

root struct Sink {
    file: io::BufferedWriter,
    workers: i64[>..] = [], data: u8[>..] = [],
    received: u8[>..] = [], offsets: i64[>..] = [], sizes: i64[>..] = [],
    blocks: i64 = 0, capacity: i64 = 0,
}

fn open(sink: Sink&, path: const u8[:], compressed: bool, threads: i64, level: i32) {
    let file = io::FileHandle { fd: 1 };
    if path.len > 0 && path != "-" {
        let opened, ok = io::open_write(path); check(ok);
        file.fd = opened.fd; opened.fd = -1;
    }
    sink.file = io::buffered_writer(file, 65536);
    if compressed {
        sink.capacity = threads * 4;
        qlimit<Block>(sink.capacity); qlimit<Encoded>(sink.capacity);
        for i in threads { sink.workers.push(thread_spawn(compress, level)); }
    }
}

fn collect(sink: Sink&) {
    sink.received.clear(); sink.offsets.clear(); sink.sizes.clear();
    sink.offsets.resize(sink.blocks, 0); sink.sizes.resize(sink.blocks, 0);
    for i in sink.blocks {
        let message = qget<Encoded>();
        sink.offsets[message.id] = sink.received.len; sink.sizes[message.id] = message.data.len;
        sink.received.append(message.data);
    }
    for i in sink.blocks {
        let bytes = sink.received[sink.offsets[i]..sink.offsets[i] + sink.sizes[i]];
        check(sink.file.write(bytes) == bytes.len);
    }
    sink.blocks = 0;
}

fn send(sink: Sink&) {
    check(qput(Block { id: sink.blocks, data: sink.data[..] }));
    sink.data.clear(); sink.blocks += 1;
    if sink.blocks == sink.capacity { sink.collect(); }
}

fn write(sink: Sink&, bytes: const u8[:]) {
    if sink.workers.len == 0 { check(sink.file.write(bytes) == bytes.len); return; }
    var at: i64 = 0;
    while at < bytes.len {
        let n = min(bytes.len - at, 65280 - sink.data.len);
        sink.data.append(bytes[at..at + n]); at += n;
        if sink.data.len == 65280 { sink.send(); }
    }
}

fn close(sink: Sink&) {
    if sink.workers.len > 0 {
        if sink.data.len > 0 { sink.send(); }
        sink.collect();
        for id in sink.workers { check(qput(Block { id: -1, data: "" })); }
        for id in sink.workers { check(thread_wait(id)); }
        sink.workers.clear();
        let eof = "\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0\x1b\0\x03\0\0\0\0\0\0\0\0\0";
        check(sink.file.write(eof) == eof.len);
    }
    check(sink.file.close());
}
