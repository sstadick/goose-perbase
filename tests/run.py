"""CLI differential tests against Rust Perbase and an independent CIGAR oracle."""
import argparse
import csv
import gzip
import io
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


def merge_tests(args, directory):
    def run(exe, data, *flags, ok=True):
        result = subprocess.run([str(exe.resolve()), "merge-adjacent", *map(str, flags)],
                                input=data, capture_output=True, timeout=60)
        assert (result.returncode == 0) == ok, (exe, flags, result.stderr)
        return result.stdout

    header = b"REF\tPOS\tEND\tDEPTH\n"
    cases = [
        (header + b"a\t0\t2\t3\na\t2\t4\t3\na\t5\t7\t3\na\t7\t8\t4\nb\t0\t1\t4\n",
         [], header + b"a\t0\t4\t3\na\t5\t7\t3\na\t7\t8\t4\nb\t0\t1\t4\n"),
        (b"a\t0\t10\t2\na\t1\t3\t2\na\t4\t6\t2", ["-n"], header + b"a\t0\t3\t2\na\t4\t6\t2\n"),
        (b"COV\tchrom\tchromStart\r\n2\ta\t0\r\n2\ta\t1\r\n", [], header + b"a\t0\t2\t2\n"),
        (b"a\t0\t\t18446744073709551615\textra\n", ["-n"], header + b"a\t0\t1\t18446744073709551615\n"),
        (b"\xef\xbb\xbf" + header + b'"a\t""b\nc"\t0\t1\t2\n', [], header + b'"a\t""b\nc"\t0\t1\t2\n'),
        (b"", [], b""), (header, [], b""), (b"\n\r\n", ["-n"], b""),
    ]
    text = io.StringIO(newline=""); writer = csv.writer(text, delimiter="\t", lineterminator="\n")
    writer.writerow(["REF", "POS", "END", "DEPTH"])
    # Cross input/output buffer boundaries, including a quoted field over 64 KiB.
    name = 'a"' + "b" * 65536 + "\nc"
    writer.writerow([name, 0, 1, 2]); writer.writerow([name, 1, 2, 2])
    for pos in range(10000): writer.writerow(["z", pos, pos + 1, pos // 7 % 3])
    oracle = io.StringIO(newline=""); writer = csv.writer(oracle, delimiter="\t", lineterminator="\n")
    writer.writerow(["REF", "POS", "END", "DEPTH"]); writer.writerow([name, 0, 2, 2])
    for pos in range(0, 10000, 7): writer.writerow(["z", pos, min(pos + 7, 10000), pos // 7 % 3])
    cases.append((text.getvalue().encode(), [], oracle.getvalue().encode()))
    for data, flags, expected in cases:
        for exe in (args.perbase, args.goose):
            assert run(exe, data, *flags) == expected, ("merge stdin", flags)
            assert run(exe, data, "-", *flags) == expected, ("merge explicit stdin", flags)
        path = directory / "merge.tsv"; path.write_bytes(data)
        assert run(args.goose, b"", path, *flags) == expected
        zipped = directory / "merge.tsv.gz"
        with pysam.BGZFile(str(zipped), "wb") as stream: stream.write(data)
        assert run(args.goose, b"", zipped, *flags) == expected
        assert run(args.perbase, b"", zipped, *flags) == expected
        packed = run(args.goose, data, *flags, "-Z", "-T", 2, "-L", 2)
        assert gzip.decompress(packed) == expected
    for data, flags in ((b"a\t0\t2\n", ["-n"]), (header + b"a\t-1\t2\t3\n", []),
                        (header + b"a\t0\t2\t18446744073709551616\n", []),
                        (header + b"a\t0\t1\n", []), (header + b"\xff\t0\t1\t2\n", [])):
        for exe in (args.goose, args.perbase): run(exe, data, *flags, ok=False)
    target = directory / "merged.tsv"
    assert run(args.goose, cases[0][0], "-o", target) == b""
    assert target.read_bytes() == cases[0][2]
    run(args.goose, b"", target, "-o", target, ok=False)
    if Path("/dev/full").exists(): run(args.goose, cases[0][0], "-o", "/dev/full", ok=False)


def cram_tests(args, directory):
    bam = directory / "cram-source.bam"
    header, records = fixture(bam)
    # Add real overlapping mates, absent SEQ and a placed unmapped read.
    for i, pos in enumerate((5, 10)):
        read = pysam.AlignedSegment()
        read.query_name = "overlap"; read.reference_id = read.next_reference_id = 0
        read.reference_start = pos; read.next_reference_start = (10, 5)[i]
        read.flag = (65, 129)[i]; read.template_length = (16, -16)[i]
        read.mapping_quality = 60; read.cigartuples = [(0, 11)]
        read.query_sequence = ("A" if i == 0 else "C") * 11; read.query_qualities = [30 + i] * 11
        records.append(read)
    for flag in (0, 4):
        read = pysam.AlignedSegment(); read.query_name = f"missing-{flag}"
        read.reference_id = 0; read.reference_start = 12; read.flag = flag
        read.cigartuples = [(0, 2), (2, 1), (3, 2), (0, 1)]; read.mapping_quality = 40
        records.append(read)
    records.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for record in records: out.write(record)
    pysam.index(str(bam))
    reference = directory / "cram-ref.fa"
    reference.write_text(">chrA\n" + "acgTN" * 19 + "A\n>chrB\n" + "TGcaN" * 6 + "C\n")
    pysam.faidx(str(reference))
    packed_reference = directory / "cram-ref.fa.gz"
    pysam.tabix_compress(str(reference), str(packed_reference), force=True)
    pysam.faidx(str(packed_reference))
    bed = directory / "cram-regions.bed"
    bed.write_text("chrB\t0\t23\nchrA\t4\t27\nchrA\t10\t15\n")

    def run(exe, command, path, flags=(), ok=True):
        result = subprocess.run([str(exe.resolve()), command, str(path), *map(str, flags)],
                                capture_output=True, timeout=60, env=dict(os.environ, RUST_LOG="error"))
        assert (result.returncode == 0) == ok, (command, path, flags, result.stderr)
        return result.stdout

    formats = [("v30", ["version=3.0"]), ("v31", ["version=3.1"]),
               ("slices", ["seqs_per_slice=3", "slices_per_container=2"]),
               ("multi", ["multi_seq_per_slice=1", "seqs_per_slice=100"]),
               ("embedded", ["embed_ref=2"]), ("noref", ["no_ref=1"])]
    for name, encoding in formats:
        cram = directory / f"{name}.cram"
        command = ["--no-PG", "-C", "-T", str(reference), "-o", str(cram)]
        for option in encoding: command += ["--output-fmt-option", option]
        pysam.view(*command, str(bam), catch_stdout=False); pysam.index(str(cram))
        # CRAM normalizes =/X CIGARs and sequence '='. Use HTSlib's decoded BAM
        # for the cross-format control, never the port's reconstructed records.
        decoded = directory / "decoded.bam"
        pysam.view("--no-PG", "-b", "-T", str(reference), "-o", str(decoded), str(cram), catch_stdout=False)
        pysam.index(str(decoded))
        for subcommand in ("base-depth", "only-depth"):
            cases = [[], ["-k", "-z", "-F", "256", "-q", "20", "-b", bed, "-M"], ["-m"]]
            cases += [["-D", "2", "-Q", "20"]] if subcommand == "base-depth" else [["-x"], ["-n", "-k"], ["-m", "--bed-format", "-z"]]
            for threads, chunk in ((1, 7), (6, 19)):
                for options in cases:
                    # Rust base-depth's bio::IndexedReader only handles plain
                    # FASTA. Keep its oracle reference uncompressed; compare
                    # Goose with the equivalent BGZF file as well.
                    oracle_flags = ["-t", threads, "-c", chunk, "-r", reference, *options]
                    flags = ["-t", threads, "-c", chunk, "-r", reference if threads == 1 else packed_reference, *options]
                    expected = run(args.perbase, subcommand, cram, oracle_flags)
                    assert run(args.goose, subcommand, decoded, flags) == expected, ("decoded BAM control", name, flags)
                    assert run(args.goose, subcommand, cram, flags) == expected, ("CRAM", name, flags)
            if subcommand == "base-depth" and name == "v30":
                flags = ["-t", 6, "-c", 19, "-r", reference]
                expected = run(args.perbase, subcommand, cram, flags)
                assert run(args.goose, subcommand, cram, [*flags, "--max-region-bytes", 1000000]) == expected
                run(args.goose, subcommand, cram, [*flags, "--max-region-bytes", 1], ok=False)
            if name in ("embedded", "noref"):
                flags = ["-t", 1, "-c", 7]
                assert run(args.goose, subcommand, cram, flags) == run(args.perbase, subcommand, cram, flags)
            else:
                run(args.goose, subcommand, cram, ["-t", 1], ok=False)
        index = Path(str(cram) + ".crai")
        index.write_bytes(gzip.decompress(index.read_bytes()))
        index.rename(cram.with_suffix(".crai"))
        flags = ["-t", 1, "-r", reference]
        assert run(args.goose, "base-depth", cram, flags) == run(args.perbase, "base-depth", decoded, flags)
        cram.with_suffix(".crai").unlink()
        run(args.goose, "base-depth", cram, flags, ok=False)
    ordinary = directory / "ordinary.fa.gz"
    ordinary.write_bytes(gzip.compress(reference.read_bytes()))
    Path(str(ordinary) + ".fai").write_bytes(Path(str(reference) + ".fai").read_bytes())
    run(args.goose, "base-depth", bam, ["-r", ordinary], ok=False)
    Path(str(packed_reference) + ".gzi").unlink()
    run(args.goose, "base-depth", bam, ["-r", packed_reference], ok=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--perbase", type=Path, required=True)
    parser.add_argument("--goose", type=Path, default=ROOT / "build/goose-perbase")
    args = parser.parse_args()
    env = dict(os.environ, RUST_LOG="error")
    with tempfile.TemporaryDirectory(prefix="goose-perbase-tests-") as t:
        directory = Path(t); bam = directory / "input.bam"
        merge_tests(args, directory)
        cram_tests(args, directory)
        def run(exe, *flags, ok=True, path=bam, threads=1, command="base-depth"):
            result = subprocess.run([str(exe.resolve()), command, str(path), "-t", str(threads), *map(str, flags)], capture_output=True, env=env, timeout=60)
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
            # Do not let HTSlib prefer a stale CSI after replacing the BAM.
            Path(str(bam) + ".csi").unlink(missing_ok=True)
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for record in sorted(records, key=lambda r: (r.reference_id, r.reference_start)): out.write(record)
            pysam.index(str(bam))

        # CLI budget reaches each worker's store; exact fit succeeds and one
        # byte less fails, independently of the read-admission depth limit.
        read = pysam.AlignedSegment()
        read.query_name = "payload"; read.reference_id = 0; read.reference_start = 2
        read.mapping_quality = 60; read.cigartuples = [(0, 1)]
        read.query_sequence = "A"; read.query_qualities = [30]; read.set_tag("ZZ", "x" * 32)
        save([read], header)
        # 7 name + 1 base + 1 quality + 4 CIGAR + (3 tag/type + 32 string + NUL).
        expected_payload = run(args.perbase)
        for threads in (1, 6):
            for limit in (49, 4294967295):
                assert run(args.goose, "--max-region-bytes", limit, threads=threads) == expected_payload
            run(args.goose, "--max-region-bytes", 48, threads=threads, ok=False)
            run(args.goose, "--max-region-bytes", 48, "-D", 0, threads=threads, ok=False)
        for limit in (-1, 0, 4294967296):
            run(args.goose, "--max-region-bytes", limit, ok=False)

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

        # only-depth: independent dense coverage oracle plus exact Rust ranges.
        def depth(exe, *flags, **kwargs):
            return run(exe, *flags, command="only-depth", **kwargs)

        for seed in (19, 47, 193):
            header, records = fixture(bam, seed)
            for fast in (False, True):
                counts = [[0] * ref["LN"] for ref in header["SQ"]]
                for r in records:
                    if r.flag & 256 or r.mapping_quality < 20: continue
                    pos = r.reference_start
                    for op, n in r.cigartuples:
                        if op in (0, 2, 3, 7, 8):
                            if op != 3 or fast:
                                for p in range(pos, pos + n): counts[r.reference_id][p] += 1
                            pos += n
                oracle = "REF\tPOS\tEND\tDEPTH\n" + "".join(
                    f'{ref["SN"]}\t{p}\t{p+1}\t{n}\n'
                    for ref, values in zip(header["SQ"], counts) for p, n in enumerate(values))
                flags = ["-n", "-k", "-z", "-F", 256, "-q", 20, *(["-x"] if fast else [])]
                assert depth(args.perbase, *flags) == oracle.encode(), ("depth oracle", seed, fast)
                assert depth(args.goose, *flags, "-c", 7, threads=6) == oracle.encode()
            # Genuine overlapping mates, including supplementary alignments and
            # long split reads. The eligibility heuristic is per record.
            for i, r in enumerate(records):
                r.query_name = f"mate-{r.reference_id}-{i // 3}"
                r.flag |= 1 | (64 if i % 2 else 128)
                r.next_reference_id = r.reference_id if i % 7 else 1 - r.reference_id
                r.next_reference_start = 10; r.template_length = -9 if i % 2 else 9
            save(records, header)
            for flags in ([], ["-x"], ["-n"], ["-n", "-k"], ["--bed-format", "-z"],
                          ["-m"], ["-m", "-n", "-F", 256, "-q", 20], ["-r", fasta],
                          ["-b", bed, "-M"], ["-B", bcf], ["-b", bed, "-B", compressed_vcf, "-M", "-m"]):
                for chunk, threads in ((7, 1), (7, 6), (1000000, 1)):
                    common = [*flags, "-c", chunk]
                    assert depth(args.goose, *common, threads=threads) == depth(args.perbase, *common, threads=threads), ("only-depth", seed, common, threads)

        # Reusing counters across long, short empty, then longer fetches must
        # clear old endpoints before shrinking and never resurrect stale depth.
        reuse_header = {"HD": {"SO": "coordinate"}, "SQ": [
            {"SN": "long", "LN": 129}, {"SN": "empty", "LN": 17}, {"SN": "longer", "LN": 150}]}
        reuse_records = []
        for tid, start, length in ((0, 96, 20), (2, 3, 4)):
            r = pysam.AlignedSegment(); r.query_name = f"reuse-{tid}"
            r.reference_id = tid; r.reference_start = start; r.cigarstring = f"{length}M"
            r.query_sequence = "A" * length; r.mapping_quality = 60; reuse_records.append(r)
        save(reuse_records, reuse_header)
        oracle = ("REF\tPOS\tEND\tDEPTH\n" + "".join(
            f'{ref["SN"]}\t{p}\t{p+1}\t{int((tid == 0 and 96 <= p < 116) or (tid == 2 and 3 <= p < 7))}\n'
            for tid, ref in enumerate(reuse_header["SQ"]) for p in range(ref["LN"]))).encode()
        for exe in (args.goose, args.perbase):
            assert depth(exe, "-n", "-k", "-z", threads=1) == oracle

        # Empty/quoted contigs still produce merged zero runs without -k.
        save([], quoted_header)
        for flags in ([], ["-n"], ["-n", "-k"], ["--bed-format", "-z"]):
            assert depth(args.goose, *flags) == depth(args.perbase, *flags)
        header, records = fixture(bam)
        records[0].query_sequence = None
        records[1].cigarstring = "1I"; records[1].query_sequence = "A"
        save(records, header)
        for flags in ([], ["-x"], ["-m"], ["-n", "-k"]):
            assert depth(args.goose, *flags) == depth(args.perbase, *flags), ("depth missing SEQ/zero span", flags)
        pysam.index("-c", str(bam)); Path(str(bam) + ".bai").unlink()
        assert depth(args.goose) == depth(args.perbase)
        pysam.index(str(bam))
        for flags in (["-m", "-x"], ["-D", 1], ["-Q", 20], ["-t", 0], ["-c", 0]):
            depth(args.goose, *flags, ok=False)
        header["SQ"][0]["LN"] = 30000; save([], header)
        plain = depth(args.perbase, "-n", "-k")
        packed = depth(args.goose, "-n", "-k", "-Z", "-T", 4, "-L", 2, "-c", 197, threads=6)
        assert gzip.decompress(packed) == plain
        assert depth(args.goose, "-n", "-k", "-o", output) == b"" and output.read_bytes() == plain
        if Path("/dev/full").exists(): depth(args.goose, "-n", "-k", "-o", "/dev/full", ok=False)
        Path(str(bam) + ".csi").unlink(missing_ok=True)
        Path(str(bam) + ".bai").unlink(); run(args.goose, ok=False)
    print("PASS: independent CIGAR/depth oracles; Rust Perbase base-depth and only-depth parity for options, workers, mates, intervals, missing SEQ, BGZF and I/O")


if __name__ == "__main__": main()
