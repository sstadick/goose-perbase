"""CLI differential tests against Rust Perbase and an independent CIGAR oracle."""
import argparse
import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile

import pysam

ROOT = Path(__file__).resolve().parents[1]
HEADER = "REF POS DEPTH A C G T N R Y S W K M INS DEL REF_SKIP FAIL COUNT_OF_MATE_RESOLUTIONS NEAR_MAX_DEPTH".split()


def fixture(path, seed=19, empty=False, depth=0):
    rng = random.Random(seed)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrA", "LN": 96}, {"SN": "chrB", "LN": 31}]}
    records = []
    shapes = [[(0, 11)], [(4, 2), (0, 3), (1, 2), (0, 4), (2, 3), (0, 2)],
              [(0, 2), (3, 5), (0, 3)], [(5, 2), (7, 3), (8, 2)], [(0, 2), (1, 4), (0, 3)]]
    for i in range(0 if empty else depth or 45):
        r = pysam.AlignedSegment()
        r.query_name = f"read-{i}"
        r.reference_id = 0 if depth else i % 2
        r.reference_start = 3 if depth else rng.randrange(15)
        r.cigartuples = [(0, 3)] if depth else shapes[i % len(shapes)]
        length = sum(n for op, n in r.cigartuples if op in (0, 1, 4, 7, 8))
        r.query_sequence = "".join(rng.choice("ACGTNMR=SVWYHKDB") for _ in range(length))
        r.query_qualities = None if i % 13 == 0 else [rng.choice([0, 19, 20, 21, 40]) for _ in range(length)]
        r.flag = 0 if depth else [0, 16, 256, 512, 1024, 2048][i % 6]
        r.mapping_quality = [0, 19, 20, 21, 60, 255][i % 6]
        records.append(r)
    records.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for record in records: out.write(record)
    pysam.index(str(path))
    return header, records


def expected(header, records, include=0, exclude=0, mapq=0, baseq=0, zeros=False, zero=False):
    counts = {}
    for r in records:
        pos, query = r.reference_start, 0
        qualities = r.query_qualities if r.query_qualities is not None else [255] * len(r.query_sequence)
        passed = (r.flag & include) == include and not r.flag & exclude and r.mapping_quality >= mapq
        for index, (op, length) in enumerate(r.cigartuples):
            if op in (0, 2, 3, 7, 8):
                for j in range(length):
                    row = counts.setdefault((r.reference_id, pos + j), [0] * 17)
                    if not passed: row[15] += 1; continue
                    if op == 3: row[14] += 1; continue
                    row[0] += 1
                    if op == 2: row[13] += 1; continue
                    base = r.query_sequence[query + j]
                    if qualities[query + j] < baseq: base = "N"
                    row[{"A": 1, "C": 2, "G": 3, "T": 4}.get(base, 5)] += 1
                    if j == length - 1 and index + 1 < len(r.cigartuples) and r.cigartuples[index + 1][0] == 1: row[12] += 1
                pos += length
            if op in (0, 1, 4, 7, 8): query += length
    lines = []
    for tid, ref in enumerate(header["SQ"]):
        for pos in range(ref["LN"]):
            if zeros or (tid, pos) in counts:
                lines.append([ref["SN"], str(pos if zero else pos + 1), *map(str, counts.get((tid, pos), [0] * 17)), "false"])
    return ("\t".join(HEADER) + "\n" + "".join("\t".join(row) + "\n" for row in lines)).encode() if lines else b""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--perbase", type=Path, required=True)
    parser.add_argument("--goose", type=Path, default=ROOT / "build/goose-perbase")
    args = parser.parse_args()
    env = dict(os.environ, RUST_LOG="error")
    with tempfile.TemporaryDirectory(prefix="goose-perbase-tests-") as t:
        directory = Path(t); bam = directory / "input.bam"
        def run(exe, *flags, ok=True, path=bam, threads=1):
            result = subprocess.run([str(exe.resolve()), "base-depth", str(path), "-t", str(threads), *map(str, flags)], capture_output=True, env=env, timeout=60)
            assert (result.returncode == 0) == ok, (exe, flags, result.stderr.decode())
            return result.stdout
        cases = [([], {}), (["-z"], {"zero": True}), (["-k"], {"zeros": True}),
                 (["-F", "3848", "-q", "20", "-Q", "20"], {"exclude": 3848, "mapq": 20, "baseq": 20}),
                 (["-f", "16", "-Q", "21"], {"include": 16, "baseq": 21}),
                 (["-q", "255", "-F", "2048", "-k", "-z"], {"mapq": 255, "exclude": 2048, "zeros": True, "zero": True})]
        for seed in (19, 47, 193):
            header, records = fixture(bam, seed)
            for flags, options in cases:
                oracle = expected(header, records, **options)
                rust = run(args.perbase, *flags)
                assert rust == oracle, ("independent oracle versus Rust", seed, flags)
                for chunk, threads in ((7, 1), (7, 6), (1000000, 1)):
                    actual = run(args.goose, *flags, "-c", chunk, threads=threads)
                    assert actual == rust, ("Goose versus Rust", seed, flags, chunk, threads)
        # CSI, basename.bai, stdout/file equality, FASTA, and empty-input semantics.
        rust = run(args.perbase)
        output = directory / "out.tsv"
        assert run(args.goose, "-o", output) == b"" and output.read_bytes() == rust
        Path(str(bam) + ".bai").rename(bam.with_suffix(".bai"))
        assert run(args.goose) == rust
        bam.with_suffix(".bai").unlink(); pysam.index("-c", str(bam))
        assert run(args.goose) == rust
        fasta = directory / "reference.fa"
        fasta.write_text(">chrA\n" + "acgTN" * 19 + "A\n>chrB\n" + "TGcaN" * 6 + "C\n")
        pysam.faidx(str(fasta))
        assert run(args.goose, "-r", fasta, "-c", 7) == run(args.perbase, "-r", fasta, "-c", 7)
        for zero in ([], ["-z"]):
            assert run(args.goose, "-r", fasta, "-k", "-c", 7, *zero, threads=6) == run(args.perbase, "-r", fasta, "-k", "-c", 7, *zero, threads=6)
        Path(str(bam) + ".csi").unlink()
        fixture(bam, empty=True)
        assert run(args.goose) == run(args.perbase) == b""
        assert run(args.goose, "-k") == run(args.perbase, "-k")
        fixture(bam, depth=199)
        assert run(args.goose, "-D", 200) == run(args.perbase, "-D", 200)
        assert b"\ttrue\n" in run(args.goose, "-D", 200)
        assert run(args.goose, "-D", 199) == run(args.perbase, "-D", 199)

        def save(records, header):
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for record in sorted(records, key=lambda r: (r.reference_id, r.reference_start)): out.write(record)
            pysam.index(str(bam))

        # Every mate strategy, including same-QNAME secondary/unpaired records,
        # must group before counting; Seqair's primary-pair links alone do not suffice.
        strategies = ["Original", "N", "IUPAC", "BaseQualMapQualFirstInPair", "BaseQualMapQualN",
                      "BaseQualMapQualIUPAC", "MapQualBaseQualFirstInPair", "MapQualBaseQualN", "MapQualBaseQualIUPAC"]
        for seed in (19, 47):
            header, records = fixture(bam, seed)
            for i, r in enumerate(records):
                r.query_name = f"group-{r.reference_id}-{i // 3}"
                r.flag |= 64 if i % 3 == 0 else 128
            save(records, header)
            for strategy in strategies:
                for filters in ([], ["-q", 20, "-Q", 21, "-F", 256]):
                    flags = ["-m", "--mate-resolution-strategy", strategy, *filters]
                    rust = run(args.perbase, *flags, "-c", 7, threads=6)
                    assert run(args.goose, *flags, "-c", 7, threads=6) == rust, ("mate strategy", seed, flags)

        # Interval union, overlap duplication, superchunk clipping and BCF RID
        # order differing from BAM order. VCF REF/END do not expand the position.
        header, records = fixture(bam)
        bed = directory / "input.bed"
        bed.write_text("chrA\t2\t16\nchrA\t0\t4\nchrA\t2\t5\nchrB\t0\t31\nchrA\t16\t21\n")
        vcf = directory / "positions.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chrB,length=31>\n##contig=<ID=chrA,length=96>\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                       "chrB\t2\t.\tAA\tA\t.\t.\t.\nchrA\t4\t.\tA\tT\t.\t.\t.\nchrA\t4\t.\tA\tG\t.\t.\t.\n")
        bcf = directory / "positions.bcf"
        with pysam.VariantFile(vcf) as src, pysam.VariantFile(bcf, "wb", header=src.header) as dst:
            for r in src: dst.write(r)
        compressed_vcf = directory / "positions.vcf.gz"
        pysam.tabix_compress(str(vcf), str(compressed_vcf))
        for selection in (["-b", bed], ["-B", vcf], ["-B", bcf], ["-B", compressed_vcf], ["-b", bed, "-B", bcf]):
            for merge in ([], ["-M"]):
                for threads in (1, 6):
                    flags = [*selection, *merge, "-k", "-z", "-c", 3, "-r", fasta]
                    rust = run(args.perbase, *flags, threads=threads)
                    assert run(args.goose, *flags, threads=threads) == rust, ("intervals", flags, threads)
        for cache in (0, 1, 10):
            assert run(args.goose, "-b", bed, "-M", "-r", fasta, "-k", "--ref-cache-size", cache, "-C", "0.000001", "-c", 3, threads=6) == run(args.perbase, "-b", bed, "-M", "-r", fasta, "-k", "-c", 3, threads=6)
        # A blank line ending exactly on a reader boundary must not hide the
        # subsequent interval. Use the same file without padding as the oracle.
        bed.write_text("chrA\t3\t8\n")
        oracle = run(args.perbase, "-b", bed, "-k")
        bed.write_bytes(b"#" + b"x" * 65533 + b"\n\nchrA\t3\t8\n")
        assert run(args.goose, "-b", bed, "-k") == oracle

        quoted_header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": 'chr"A', "LN": 3}]}
        save([], quoted_header)
        assert run(args.goose, "-k") == run(args.perbase, "-k")

        # Read-admission boundaries: ending exactly at a new start, new-start
        # bypass, excluded reads occupying capacity, and chunk-local replay.
        header, records = fixture(bam, depth=5)
        for i, r in enumerate(records):
            r.reference_start = [0, 0, 1, 1, 2][i]
            r.cigarstring = ["1M", "3M", "2M", "1M", "1M"][i]
            r.query_sequence = "ACGTA"[i] * r.infer_query_length(); r.query_qualities = [30] * len(r.query_sequence)
            r.flag = 256 if i == 0 else 0
        save(records, header)
        for depth in (0, 1, 2, 3, 100000):
            for chunk in (1, 3, 100):
                flags = ["-D", depth, "-c", chunk, "-F", 256, "-k"]
                assert run(args.goose, *flags, threads=6) == run(args.perbase, *flags, threads=6), ("admission", depth, chunk)
        records[0].cigarstring = "1M1I1M1D1N1M"; records[0].query_sequence = None; records[0].flag = 0
        save(records[:1], header)
        for flags in ([], ["-Q", 0], ["-Q", 20], ["-m"], ["-k", "-r", fasta]):
            assert run(args.goose, *flags) == run(args.perbase, *flags), ("missing SEQ", flags)

        # Zero-span reads create no column but can still occupy admission slots.
        # Dropping the 1I before admission incorrectly keeps the final 1M at D=1.
        for i, r in enumerate(records[:3]):
            r.reference_start = [0, 1, 1][i]; r.cigarstring = ["1M", "1I", "1M"][i]
            r.query_sequence = "A"; r.query_qualities = [30]; r.flag = 0
        save(records[:3], header)
        for depth in (0, 1, 2):
            for chunk in (1, 100):
                flags = ["-D", depth, "-c", chunk, "-k", "-z"]
                assert run(args.goose, *flags) == run(args.perbase, *flags), ("zero-span admission", depth, chunk)

        # Enough output for many BGZF blocks/batches; decompress externally and
        # assert ordering/content independently of the compression thread count.
        header["SQ"][0]["LN"] = 30000
        save([], header)
        plain = run(args.perbase, "-k")
        for threads, level in ((1, 0), (2, 2), (4, 9)):
            packed = run(args.goose, "-k", "-Z", "-T", threads, "-L", level, "-c", 197, threads=6)
            assert gzip.decompress(packed) == plain
            assert packed[-28:] == bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
        for flags in (["-t", 0], ["-c", 0], ["-q", 256], ["-F", -1], ["--mate-resolution-strategy", "invalid"], ["-T", 0], ["-L", 10], ["-C", "nan"]):
            run(args.goose, *flags, ok=False)
        if Path("/dev/full").exists():
            run(args.goose, "-k", "-o", "/dev/full", ok=False)
            run(args.goose, "-k", "-Z", "-o", "/dev/full", ok=False)
        Path(str(bam) + ".bai").unlink(); run(args.goose, ok=False)
    print("PASS: independent CIGAR oracle; Rust Perbase parity for threads, all mate strategies, BED/VCF/BCF, depth caps, FASTA/zeros, missing SEQ, BGZF and I/O")


if __name__ == "__main__": main()
