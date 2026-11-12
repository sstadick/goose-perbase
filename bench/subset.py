"""Extract an indexed BAM region for repeatable end-to-end benchmarks; keep original coordinates/header."""
import argparse
from pathlib import Path
import pysam

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input", type=Path)
p.add_argument("output", type=Path)
p.add_argument("contig")
p.add_argument("begin", type=int)
p.add_argument("end", type=int)
args = p.parse_args()
if args.output.exists(): p.error("output already exists")
if args.begin < 0 or args.end <= args.begin: p.error("expected nonempty zero-based half-open region")
with pysam.AlignmentFile(args.input, "rb") as source:
    with pysam.AlignmentFile(args.output, "wb", header=source.header) as output:
        count = 0
        for record in source.fetch(args.contig, args.begin, args.end):
            output.write(record); count += 1
pysam.index(str(args.output))
print(f"{count} records, {args.output.stat().st_size} compressed bytes")
