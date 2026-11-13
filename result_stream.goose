// Transfer region output in bounded pieces, preserving the main thread's
// region order. Goose's typed queues cannot wait for a particular region;
// the native mailboxes provide that wait without polling. Jobs still travel
// through ordinary Goose queues.
namespace result_stream;
import std;

// Each mailbox holds at most capacity bytes, with a minimum capacity of 64 KiB.
// Workers can produce ahead of the writer until their mailbox fills. Storage
// is allocated on first use and reused for later regions. The main thread owns
// the handle and closes it only after all workers have joined. Workers share
// this opaque handle, never a Goose reference or a BAM reader. A slot belongs
// to one outstanding region; reuse it only after receive() returns zero for it.
extern fn open(slots: i64, capacity: i64) -> u64 = "perbase_stream_open";
extern fn close(handle: u64) = "perbase_stream_close";
extern fn send(handle: u64, slot: i64, bytes: const u8[:]) = "perbase_stream_send";

// Copy the next piece into a 64 KiB destination and return its byte count.
// Zero marks the end of a region, including a region with no output rows.
extern fn receive(handle: u64, slot: i64, bytes: u8[:]) -> i64 = "perbase_stream_receive";

// Call after appending each row. Keep at most one partial piece plus the
// current row in text. Large rows are split across messages, not truncated.
// Sending blocks only when the region's mailbox has no room for another piece,
// so later regions cannot accumulate unlimited out-of-order output.
fn flush(handle: u64, slot: i64, text: u8[>..]&, finish: bool = false) {
    if text.len < 65536 && !finish {
        return;
    }
    var at: i64 = 0;
    while at < text.len {
        let end = min(at + 65536, text.len);
        send(handle, slot, text[at..end]);
        at = end;
    }
    text.clear();
    if finish {
        send(handle, slot, "");
    }
}
