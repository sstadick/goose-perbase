# goose-perbase

An experimental Goose port of [Perbase](https://github.com/sstadick/perbase),
implementing **`base-depth`: indexed BAM → TSV or BGZF-compressed TSV**. It uses
[goose-seqair](https://github.com/sstadick/goose-seqair)'s indexed region slabs
and raw column iterator, rather than an HTSlib pileup wrapper. The application
and aggregation are Goose; BGZF decompression and base unpacking use the native
kernels supplied by its dependencies. No Rust or HTSlib runtime is needed.

**[Benchmarks and reproduction](bench/README.md)** · **[Aggregation API](API.md)**

## Build and run

Requires Goose, Python 3.11+, Clang, and CMake 3.20+. Tested with Goose
[`74c641c5`](https://github.com/aardappel/goose/commit/74c641c5ac6d762e57e3177678fe9dea414d042a)
and Clang 18.1.3 on Linux x86-64. Runtime imports pin released dependencies,
including goose-seqair v0.1.0 (zero-span retention for depth-cap parity) and
goose-bio v0.3.1 (indexed FASTA), with goose-io v0.3.0 for file seeking.
No sibling library checkout is needed.

```sh
python3 scripts/build.py
build/goose-perbase base-depth reads.bam -o depth.tsv
build/goose-perbase base-depth reads.bam -F 3848 -q 20 -Q 20 -o filtered.tsv
build/goose-perbase base-depth reads.bam -r reference.fa -o with-reference.tsv
build/goose-perbase base-depth reads.bam -t 8 -m -b targets.bed -o targets.tsv
build/goose-perbase base-depth reads.bam -t 8 -Z -T 4 -L 2 -o depth.tsv.gz
build/goose-perbase base-depth --help
```

The build compiles generated C, native adapters, and libdeflate using
**Clang `-O3 -march=native`**. Goose's `--standalone -O2` is its source-level
inlining setting, not the native compiler's optimization level. Native binaries
target the build machine's CPU; rebuild on a different machine.

The BAM must be coordinate-sorted and indexed (`reads.bam.bai`,
`reads.bam.csi`, or `reads.bai`). A reference is optional and must be
uncompressed FASTA with a `.fai` index. Output defaults to stdout; `-o -` also
selects stdout. Compression requires `-Z`, not a filename suffix. Errors exit
nonzero; a failure after processing begins may leave partial output, which must
be discarded.

## Compatibility and current limits

The compatibility target is **Perbase v1.4.1**, revision
[`213e774`](https://github.com/sstadick/perbase/commit/213e77439b2c216662528c5eebfff2d18492d025).
All its `base-depth` processing option names are accepted within the indexed
BAM input scope. Other Perbase subcommands, SAM/CRAM/stdin input, and compressed
FASTA are not implemented.

| Options | Behavior/defaults |
| --- | --- |
| `-o`, `--output` | File or stdout (`-`) |
| `-t`, `--threads` | Available CPUs; clamp request to available CPUs, then use `max(threads - 2, 1)` pileup workers |
| `-c`, `--chunksize` | 1,000,000 bases |
| `-C`, `--channel-size-modifier` | 0.15; approximate output-buffer GB per worker |
| `-f`, `--include-flags`; `-F`, `--exclude-flags` | Both 0 |
| `-q`, `--min-mapq`; `-Q`, `--min-base-quality-score` | No quality filtering by default |
| `-m`, `--mate-fix`; `--mate-resolution-strategy` | Off; strategy `Original` |
| `-b`, `--bed-file`; `-B`, `--bcf-file` | Union of BED intervals and VCF/BCF positions |
| `-M`, `--skip-merging-intervals` | Off; preserve duplicate rows when enabled |
| `-r`, `--ref-fasta`; `--ref-cache-size` | Optional reference; 10 cached windows per worker |
| `-k`, `--keep-zeros`; `-z`, `--zero-base` | Both off |
| `-D`, `--max-depth` | 100,000; admission before filtering, not a per-position clamp |
| `-Z`, `--bgzip`; `-T`, `--compression-threads`; `-L`, `--compression-level` | Off; 4 compression workers; level 2 (0–9 accepted) |

Mate strategies are case-insensitive: `Original`, `N`, `IUPAC`,
`BaseQualMapQualFirstInPair`, `BaseQualMapQualN`, `BaseQualMapQualIUPAC`,
`MapQualBaseQualFirstInPair`, `MapQualBaseQualN`, and `MapQualBaseQualIUPAC`.
Grouping uses query names, including secondary and unpaired alignments, rather
than only primary-pair links. `-M` is exclusively the interval-merging option.

TSV uses the corrected `COUNT_OF_MATE_RESOLUTIONS` spelling. Excluded reads
contribute to `FAIL` rather than disappearing during BAM loading. Deletions
count toward depth; reference skips do not. Insertions count events, not inserted bases. Low-quality
bases, missing SEQ, and input IUPAC ambiguity codes count as N. Reference bases
are included for synthesized zero rows too. An empty BAM produces no
header unless zero rows are requested. Without mate fixing, the mate-resolution
counter and ambiguity-resolution columns remain zero.

Intervals sort by BAM reference order, start, then end. Touching intervals merge
by default; VCF/BCF records contribute POS only, ignoring REF length/INFO END.
BED/VCF text and BCF 2.2 support plain or BGZF input, not ordinary gzip streams.
Without merging, output preserves interval concatenation within each Perbase
superchunk (`chunksize × workers`), so positions can repeat or move backward.

The depth cap applies to whole reads before read filters and mate resolution.
The first read at a later start can bypass the cap; depth may therefore exceed
`-D`. Admission restarts per fetch, so changing chunk boundaries can change
capped results. `-D 0` is not unlimited.

Zero-reference-span alignments remain available for admission but produce no
columns. Perbase's mate comparator is not a total ordering for some
multi-alignment groups; matching small-group insertion order does not guarantee identical results for arbitrary
large groups. Missing-SEQ overlapping mates are handled safely here, whereas
some upstream strategies index absent sequence/qualities.

Goose-seqair's default 256 MiB variable-payload limit applies to each loaded
region; metadata and reusable allocation capacity add memory beyond that.
Reduce `-c` for unusually dense regions. Restricted intervals may span an entire
superchunk rather than one chunk. BAM/index consistency and sorting in
skipped regions are not exhaustively validated.

## Parallel execution and memory

Persistent Goose `thread_fn` workers each open and reuse their own BAM reader,
parsed index, region store, reference handle/cache, and output buffers. No file
descriptor or borrowed view crosses a thread boundary. Jobs and results contain
owned values copied through typed queues. Ordered, bounded batches keep output
deterministic without Rayon or a shared reader lock. Compression uses a separate
persistent worker pool and writes completed BGZF blocks in submission order.

Workers serialize each region into TSV; the coordinator collects a bounded batch
and writes it in interval order. `-C` limits concurrent jobs using an estimated
serialized-row size, not Rust's position-struct size, and is **not a total-memory
limit**. BAM slabs, copies in queues, reference buffers and reused capacities
add memory. Compression queues hold at most four blocks per compression worker.
`--ref-cache-size` bounds reference windows per worker, unlike upstream's shared
full-contig cache; the tuning knob has different memory/performance semantics.

## Validate

The differential suite requires a separately built Rust Perbase v1.4.1 binary;
see the benchmark instructions for its exact revision and native build flags.
Python/pysam creates fixtures and supplies an independent CIGAR-counting oracle;
it is not used by the executable or timed benchmark commands.

```sh
python3 -m venv .ci/venv
.ci/venv/bin/pip install -r tests/requirements.txt
.ci/venv/bin/python tests/run.py --perbase /path/to/perbase
python3 scripts/build.py --check
.ci/venv/bin/python tests/run.py --perbase /path/to/perbase --goose build/goose-perbase-check
```

The checked build enables Goose bounds checks and ASan/UBSan. Cases cover worker
counts, all mate strategies, BED/VCF/BCF union/duplicates, superchunk boundaries,
CIGAR operations, flags/qualities, missing SEQ, BAI/CSI, FASTA/zero rows, depth
admission, quoted names, reader boundaries, empty input, and output failures.
BGZF output is decompressed independently. Real-data benchmarks additionally
compare every timed TSV byte-for-byte by SHA-256.

API documentation follows Gooseverse's goose-doc conventions. With goose-doc
checked out next door (or `GOOSE_DOC` pointing to its `main.goose`):

```sh
CC=clang CXX=clang++ GOOSE_CC=clang scripts/docs.sh --check
CC=clang CXX=clang++ GOOSE_CC=clang scripts/doctest.sh
```

MIT licensed; counting semantics and CLI names follow the original Perbase.
