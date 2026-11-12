# Base-depth benchmarks

## Current comparison: Perbase v1.4.1

The expanded CLI is compared against Rust Perbase
[`213e774`](https://github.com/sstadick/perbase/commit/213e77439b2c216662528c5eebfff2d18492d025)
(v1.4.1), which **already reuses indexed readers per worker**. Do not use the
historical v1.1.0 speedups below to describe this implementation.

Median wall time; memory is the highest measured RSS. All runs produced matching
TSV bytes (after decompression for BGZF).

| Paper chr20, 10 Mb | Goose seconds | Perbase seconds | Speedup | Goose MiB | Perbase MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| [Plain TSV, `-t 1`](results/current-paper-t1.json) | 5.877 | 8.476 | **1.44×** | 415.4 | 474.3 |
| [Plain TSV, `-t 8`](results/current-paper-t8.json) | 3.073 | 5.430 | **1.77×** | 1937.7 | 2186.0 |
| [Mate fixing, `-t 8 -m`](results/current-paper-mates-t8.json) | 3.645 | 18.849 | **5.17×** | 1936.1 | 2193.1 |
| [BGZF, `-t 8 -Z -T 4 -L 2`](results/current-paper-bgzf-t8.json) | 3.047 | 4.913 | **1.61×** | 1918.7 | 2193.4 |

Six Goose pileup workers reduce plain-TSV runtime by **1.91×** versus one worker,
not sixfold. The default runs emit 569,956,346 TSV bytes each. These are complete
application timings, including serialization, copied results and ordered output;
they do not establish which stage limits scaling. Compared with the historical
Goose baseline, region buffering also increases single-worker memory and runtime.

Both programs process the paper chr20 10 Mb subset described below, using
`-c 1000000` and identical CPU affinity `0,1,2,3,4,5,6,7` (eight physical cores
on this Ryzen 9 3950X). `-t 1` means one pileup worker and `-t 8` means six.
Plain TSV comparisons use seven measured runs; mate-fixing and BGZF controls
use three, after one warmup each. Orders alternate; every output is checked by
SHA-256 outside timing. BGZF comparisons hash the independently decompressed
TSV, not compressor-dependent bytes. Output is written to real files under
`/srv`, including flush/close but excluding `fsync`.

Build settings remain Clang `-O3 -march=native` for generated C and C libraries,
and Rust 1.98.1 with `-C opt-level=3 -C target-cpu=native -C linker=clang`.
The Goose compiler is unchanged. The application uses the sibling Seqair
checkout at `35266830e97e` plus its unpublished zero-span retention change;
JSON reports include the application and modified dependency source hashes,
binary hashes, commands, input hash, output hash, samples and peak RSS.

Reproduce after building both applications as described below, substituting
the v1.4.1 revision and a fresh Rust target directory:

```sh
export TMPDIR=/srv/goose-perbase
export PERBASE=/srv/goose-perbase/rust-current-target/release/perbase
for threads in 1 8; do
  python3 bench/run.py "$TMPDIR/paper-chr20-10mb.bam" \
    --perbase "$PERBASE" --cpus 0,1,2,3,4,5,6,7 --threads "$threads" \
    --rounds 7 --output "$TMPDIR/current-paper-t${threads}.json"
done
python3 bench/run.py "$TMPDIR/paper-chr20-10mb.bam" \
  --perbase "$PERBASE" --cpus 0,1,2,3,4,5,6,7 --threads 8 \
  --rounds 3 --output "$TMPDIR/current-paper-mates-t8.json" --extra -m
python3 bench/run.py "$TMPDIR/paper-chr20-10mb.bam" \
  --perbase "$PERBASE" --cpus 0,1,2,3,4,5,6,7 --threads 8 \
  --rounds 3 --output "$TMPDIR/current-paper-bgzf-t8.json" --extra -Z -T 4 -L 2
```

The full 35 GiB whole-genome BAM has not been benchmarked. Current Goose workers
buffer serialized regions and send owned results in bounded, ordered batches;
they do not stream individual rows directly to the output as the old baseline
did. Parallelism trades additional slab/output-buffer memory for throughput.

## Historical single-worker baseline: Perbase v1.1.0

The following measurements belong to the earlier, restricted implementation
and the older Rust baseline. Commands and implementation descriptions in this
section describe that historical comparison, not the current application.

These are end-to-end **indexed BAM → uncompressed TSV** measurements, including
process startup, index/header loading, decompression, pileups, aggregation, and
output. They do not isolate the pileup engine. Both programs use `-t 1 -c 1000000`,
default filters, no mate fixing, no reference, and no zero filling.

`run.py` warms each program once, alternates their order across seven measured
runs, and verifies every completed TSV against the first output using SHA-256
outside the timed interval. It writes real files, not `/dev/null`. Times include
buffer flushing and close, but not an `fsync`: these are warm-cache filesystem
measurements, not durable-storage throughput. `/usr/bin/time` supplies peak RSS.
The runner fails on nonzero exits or differing output.

### Historical results

Median wall time across seven measured runs; memory is the highest RSS across
those runs. All four inputs produced matching TSV hashes on every run.

| Input | Goose seconds | Perbase seconds | Speedup | Goose MiB | Perbase MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| [Paper chr20, 10 Mb; 2.24M records](results/paper-chr20-10mb.json) | 4.374 | 21.103 | **4.82×** | 141.2 | 475.1 |
| [Paper tiny; 217K records](results/paper-tiny.json) | 0.520 | 10.234 | 19.69× | 161.2 | 343.3 |
| [High-depth subset; 101K records](results/spike.json) | 0.191 | 0.764 | 4.00× | 55.5 | 186.5 |
| [Seqair test BAM; 29K records](results/initial-small.json) | 0.039 | 0.185 | 4.74× | 8.0 | 185.9 |

The largest comparison emitted 9,999,137 data rows and 569,956,345 identical TSV
bytes per run. A separate Goose bounds-checked ASan/UBSan run also produced the
same hash without sanitizer errors. These are single-worker subset results,
**not whole-genome or multi-core results**.
The tiny fixture retains a large reference dictionary despite sparse coverage;
its 19.69× result should not be generalized (see Interpretation below).

### Historical build provenance

Measurements made on 2026-11-12: AMD Ryzen 9 3950X, Linux x86-64, Clang 18.1.3,
with both processes and their child threads pinned to logical CPU 4. Native
code and C dependencies use **Clang `-O3 -march=native`**; Rust uses optimization
level 3, `target-cpu=native`, and Clang as linker. No GCC builds are compared.

- Goose compiler:
  [`74c641c5ac6d762e57e3177678fe9dea414d042a`](https://github.com/aardappel/goose/commit/74c641c5ac6d762e57e3177678fe9dea414d042a).
- Goose-seqair:
  [`35266830e97e0c7f06cdc751f66ae5ed95286232`](https://github.com/sstadick/goose-seqair/commit/35266830e97e0c7f06cdc751f66ae5ed95286232).
- Rust Perbase 1.1.0:
  [`c761eb9b4ce53eb6e7110dcb882ee5627b97b0c6`](https://github.com/sstadick/perbase/commit/c761eb9b4ce53eb6e7110dcb882ee5627b97b0c6),
  built with Rust 1.98.1 and its committed Cargo.lock.

From goose-perbase:

```sh
python3 scripts/build.py
```

From a Perbase checkout at the revision above, use a fresh target directory to
avoid stale native objects from different compiler settings:

```sh
CC=clang CXX=clang++ \
CFLAGS='-O3 -march=native' CXXFLAGS='-O3 -march=native' \
RUSTFLAGS='-C opt-level=3 -C target-cpu=native -C linker=clang' \
CARGO_TARGET_DIR=/srv/goose-perbase/rust-target \
cargo +1.98.1 build --release --locked
```

`scripts/build.py` emits standalone C with Goose's `-O2` frontend inlining
setting, then compiles it and native sources with Clang `-O3 -march=native`.
It also rebuilds libdeflate with those native settings. `run.py` records binary
hashes but does not establish how an externally supplied binary was built;
reproduce the builds above before comparing timings.

### Inputs and historical reproduction

Run these commands from goose-perbase. Python/pysam is only needed to extract
subsets; see `tests/requirements.txt`. Point `PYTHON` to that environment, and
`PERBASE` to the freshly built Rust binary:

```sh
export TMPDIR=/srv/goose-perbase
mkdir -p "$TMPDIR"
export PERBASE=/srv/goose-perbase/rust-target/release/perbase
export PYTHON=.ci/venv/bin/python
```

The Perbase paper uses 1000 Genomes sample HG00157, approximately 30× whole-genome
coverage. Its [data-fetch script](https://github.com/sstadick/perbase/blob/c761eb9b4ce53eb6e7110dcb882ee5627b97b0c6/paper/fetch_data.sh)
documents the original CRAM/reference sources and BAM conversion. We used the
existing indexed BAMs from that checkout; benchmark outputs never overwrite
the paper's saved outputs.

#### Paper tiny fixture

`HG00157.tiny.bam` contains 217,454 records (216,860 mapped), 18,184,667 compressed
bytes, and the full 3,366-reference dictionary. It is sparse relative to that
dictionary: mapped records occur only on chr1.

```sh
python3 bench/run.py ../perbase/paper/data/HG00157.tiny.bam \
  --perbase "$PERBASE" --cpu 4 --output "$TMPDIR/paper-tiny.json"
```

#### Paper chr20, 10 Mb subset

Extract records overlapping zero-based **chr20:[10,000,000,20,000,000)** from
`HG00157.final.bam`, retaining original headers, coordinates, and complete reads.
Reads crossing the boundaries are not trimmed, so output can extend slightly
beyond the interval. Both programs process the same extracted BAM, not a
region-limited run on different inputs. Extraction/indexing is outside timing.

```sh
"$PYTHON" bench/subset.py ../perbase/paper/data/HG00157.final.bam \
  "$TMPDIR/paper-chr20-10mb.bam" chr20 10000000 20000000
python3 bench/run.py "$TMPDIR/paper-chr20-10mb.bam" \
  --perbase "$PERBASE" --cpu 4 --output "$TMPDIR/paper-chr20-10mb.json"
```

The extracted BAM contains 2,241,960 records (2,239,035 mapped), 109,209,558
compressed bytes, and the original 3,366-reference dictionary. Its source BAM
is 36,861,040,449 bytes. The full whole-genome BAM has **not** been benchmarked.

#### Additional controls

- Seqair's `test.bam`: 29,358 records, 1,581,777 bytes. Pass it directly to
  `bench/run.py` with the same options.
- Local `spike1.bam`, chr1:[155,850,000,155,950,000): extract with
  `bench/subset.py` just as above. This high-depth subset contains 101,277
  records, 3,356,997 bytes. This local input is not distributed with the project.

The JSON reports retain input, binary, and output hashes; exact compressed bytes
may differ if a subset is regenerated with another pysam/HTSlib version.

### Historical interpretation

This measured the earlier applications, not a drop-in feature-complete Perbase
replacement. See the [scope and limitations](../README.md#compatibility-and-current-limits).
There is no multi-core comparison or extrapolation to whole-genome runtime.

Perbase opens a new indexed reader and collects a vector of positions for each
region. Goose reuses the reader/index/slabs and streams TSV rows through reusable
buffers. Sparse data retaining a whole-genome header can emphasize setup costs;
a large tiny-fixture speedup must not be presented as an isolated pileup-kernel
speedup. Neither implementation has been modified specifically for these data.
