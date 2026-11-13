#include <assert.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Goose's byte-slice ABI is packed. No pointer into a Goose buffer is retained. */
#pragma pack(push, 1)
typedef struct {
    uint8_t *data;
    int64_t len;
} perbase_bytes;
#pragma pack(pop)

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t readable;
    pthread_cond_t writable;
    size_t head;
    size_t used;
    int finished;
    uint8_t *data;
} mailbox;

typedef struct {
    int64_t count;
    size_t capacity;
    mailbox slots[];
} stream;

static void checked(int error) {
    if (error != 0) {
        abort();
    }
}

uint64_t perbase_stream_open(int64_t count, int64_t capacity) {
    if (count <= 0 || (uint64_t)count > (SIZE_MAX - sizeof(stream)) / sizeof(mailbox) ||
        capacity < 65536 || (uint64_t)capacity > SIZE_MAX) {
        return 0;
    }
    stream *s = calloc(1, sizeof(stream) + (size_t)count * sizeof(mailbox));
    if (!s) {
        return 0;
    }
    s->count = count;
    s->capacity = (size_t)capacity;
    for (int64_t i = 0; i < count; ++i) {
        checked(pthread_mutex_init(&s->slots[i].mutex, NULL));
        checked(pthread_cond_init(&s->slots[i].readable, NULL));
        checked(pthread_cond_init(&s->slots[i].writable, NULL));
    }
    return (uint64_t)(uintptr_t)s;
}

void perbase_stream_send(uint64_t handle, int64_t index, perbase_bytes bytes) {
    stream *s = (stream *)(uintptr_t)handle;
    assert(index >= 0 && index < s->count && bytes.len >= 0 && bytes.len <= 65536);
    mailbox *slot = &s->slots[index];
    checked(pthread_mutex_lock(&slot->mutex));
    assert(!slot->finished);
    if (bytes.len) {
        while (s->capacity - slot->used < (size_t)bytes.len) {
            checked(pthread_cond_wait(&slot->writable, &slot->mutex));
        }
        /* Empty regions need no allocation. Reuse storage across later jobs. */
        if (!slot->data) {
            slot->data = malloc(s->capacity);
            if (!slot->data) {
                abort();
            }
        }
        size_t tail = (slot->head + slot->used) % s->capacity;
        size_t first = s->capacity - tail;
        if (first > (size_t)bytes.len) {
            first = (size_t)bytes.len;
        }
        memcpy(slot->data + tail, bytes.data, first);
        memcpy(slot->data, bytes.data + first, (size_t)bytes.len - first);
        slot->used += (size_t)bytes.len;
    } else {
        /* The reader returns the end marker only after draining these bytes. */
        slot->finished = 1;
    }
    checked(pthread_cond_signal(&slot->readable));
    checked(pthread_mutex_unlock(&slot->mutex));
}

int64_t perbase_stream_receive(uint64_t handle, int64_t index, perbase_bytes bytes) {
    stream *s = (stream *)(uintptr_t)handle;
    assert(index >= 0 && index < s->count && bytes.len == 65536);
    mailbox *slot = &s->slots[index];
    checked(pthread_mutex_lock(&slot->mutex));
    while (!slot->used && !slot->finished) {
        checked(pthread_cond_wait(&slot->readable, &slot->mutex));
    }
    size_t length = slot->used < 65536 ? slot->used : 65536;
    if (length) {
        size_t first = s->capacity - slot->head;
        if (first > length) {
            first = length;
        }
        memcpy(bytes.data, slot->data + slot->head, first);
        memcpy(bytes.data + first, slot->data, length - first);
        slot->head = (slot->head + length) % s->capacity;
        slot->used -= length;
    } else {
        /* The coordinator may now assign another job to this empty slot. */
        slot->head = 0;
        slot->finished = 0;
    }
    checked(pthread_cond_signal(&slot->writable));
    checked(pthread_mutex_unlock(&slot->mutex));
    return (int64_t)length;
}

/* All users have joined before destruction. Fatal worker/I/O errors terminate
 * the process, as in the Goose pipeline; there is no detached producer. */
void perbase_stream_close(uint64_t handle) {
    stream *s = (stream *)(uintptr_t)handle;
    for (int64_t i = 0; i < s->count; ++i) {
        checked(pthread_cond_destroy(&s->slots[i].readable));
        checked(pthread_cond_destroy(&s->slots[i].writable));
        checked(pthread_mutex_destroy(&s->slots[i].mutex));
        free(s->slots[i].data);
    }
    free(s);
}
