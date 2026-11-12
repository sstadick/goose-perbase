"""Alternate whole-process Perbase runs, validating every TSV before reporting speed."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def digest(path, compressed=False):
    with (gzip.open(path, "rb") if compressed else path.open("rb")) as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bam", type=Path, help="Indexed alignment input, or TSV for merge-adjacent")
    p.add_argument("--perbase", type=Path, required=True)
    p.add_argument("--goose", type=Path, default=ROOT / "build/goose-perbase")
    p.add_argument("--rounds", type=int, default=7)
    p.add_argument("--cpu", type=int)
    p.add_argument("--cpus", help="Comma-separated logical CPUs for a parallel comparison")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--command", choices=("base-depth", "only-depth", "merge-adjacent"), default="base-depth")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--rust-reference", type=Path, help="Equivalent plain FASTA for Rust's base-depth reader; overrides -r/--ref-fasta in --extra")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="Common additional subcommand arguments")
    args = p.parse_args()
    if args.rounds < 1: p.error("rounds must be positive")
    if args.threads < 1: p.error("threads must be positive")
    if args.cpu is not None and args.cpus: p.error("choose --cpu or --cpus")
    if args.cpu is not None: os.sched_setaffinity(0, {args.cpu})
    if args.cpus: os.sched_setaffinity(0, set(map(int, args.cpus.split(','))))
    compressed = "-Z" in args.extra or "--bgzip" in args.extra
    report = {"platform": platform.platform(), "affinity": sorted(os.sched_getaffinity(0)),
              "clang": subprocess.check_output(["clang", "--version"], text=True).splitlines()[0],
              "input_sha256": digest(args.bam), "input_bytes": args.bam.stat().st_size,
              "rounds": args.rounds, "command": args.command, "threads_flag": None if args.command == "merge-adjacent" else args.threads, "extra_flags": args.extra, "native_flags": "-O3 -march=native",
              "note": "Build both binaries with documented flags first. This runner does not establish external binary provenance.",
              "results": {}}
    report["source_sha256"] = {path.name: digest(path) for path in sorted(ROOT.glob("*.goose"))}
    packages = Path(os.environ.get("GOOSE_HOME", Path.home() / ".goose")) / "pkg/github.com/sstadick"
    report["seqair_source_sha256"] = {name: digest(packages / "goose-seqair/v0.2.2" / name)
                                      for name in ("pileup.goose", "cram.goose", "cram_index.goose", "src/cram.c")}
    report["bio_source_sha256"] = {"faidx.goose": digest(packages / "goose-bio/v0.4.1/faidx.goose")}
    rust_extra = args.extra.copy()
    if args.rust_reference:
        references = [i for i, flag in enumerate(rust_extra) if flag in ("-r", "--ref-fasta")]
        if len(references) != 1 or references[0] + 1 >= len(rust_extra):
            p.error("--rust-reference requires exactly one -r/--ref-fasta PATH in --extra")
        rust_extra[references[0] + 1] = str(args.rust_reference.resolve())
        report["rust_reference_override"] = str(args.rust_reference.resolve())
    programs = {"goose": args.goose.resolve(), "rust": args.perbase.resolve()}
    samples = {name: [] for name in programs}; rss = {name: [] for name in programs}
    expected = None
    with tempfile.TemporaryDirectory(prefix="perbase-bench-") as t:
        directory = Path(t)
        for round_ in range(args.rounds + 1):
            for name in programs if round_ % 2 == 0 else reversed(programs):
                output = directory / f"{name}.tsv"; memory = directory / f"{name}.rss"
                depth_flags = [] if args.command == "merge-adjacent" else ["-t", str(args.threads), "-c", "1000000"]
                extra = rust_extra if name == "rust" else args.extra
                command = [str(programs[name]), args.command, str(args.bam.resolve()), *depth_flags, "-o", str(output), *extra]
                report.setdefault("commands", {})[name] = command
                start = time.perf_counter_ns()
                subprocess.run(["/usr/bin/time", "-f", "%M", "-o", str(memory), *command],
                               env=dict(os.environ, RUST_LOG="error"), check=True, stdout=subprocess.DEVNULL)
                elapsed = (time.perf_counter_ns() - start) / 1e9
                result = digest(output, compressed)
                if expected is None: expected = result
                assert result == expected, (name, "TSV differs", expected, result)
                if round_ > 0: samples[name].append(elapsed); rss[name].append(int(memory.read_text()))
                report.setdefault("encoded_output_bytes", {})[name] = output.stat().st_size
                print(f"{'warmup' if round_ == 0 else round_} {name}: {elapsed:.6f}s; output verified", flush=True)
        report["output_sha256"] = expected
        for name, exe in programs.items():
            report["results"][name] = {"binary_sha256": digest(exe), "samples_seconds": samples[name],
                "median_seconds": statistics.median(samples[name]), "peak_rss_kib": max(rss[name])}
        report["speedup"] = statistics.median(samples["rust"]) / statistics.median(samples["goose"])
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
