"""Build all timed native code with Clang -O3 -march=native, including libdeflate."""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Debug range checks plus ASan/UBSan; writes goose-perbase-check")
    args = parser.parse_args()
    build = ROOT / "build"; build.mkdir(exist_ok=True)
    source = build / ("check.c" if args.check else "main.c")
    env = dict(os.environ, CC="clang", CXX="clang++", GOOSE_CC="clang")
    subprocess.run([env.get("GOOSE", "goose"), "--standalone", "-O2", *(["--debug"] if args.check else []),
                    "-o", str(source), "main.goose"], cwd=ROOT, env=env, check=True)
    headers = [(source.parent / p).resolve() for p in re.findall(r'^#include "([^"]+)"', source.read_text(), re.M)]
    def package(header): return next(p.parent.parent for p in headers if p.name == header)
    codec = package("deflate.h")
    subprocess.run(["cmake", "-S", str(codec / "vendor/libdeflate"), "-B", str(build / "codec"),
                    f"-DCMAKE_C_COMPILER={shutil.which('clang')}", "-DCMAKE_BUILD_TYPE=Release",
                    "-DCMAKE_C_FLAGS_RELEASE=-O3 -march=native -DNDEBUG", "-DLIBDEFLATE_BUILD_SHARED_LIB=OFF",
                    "-DLIBDEFLATE_BUILD_GZIP=OFF", "-DLIBDEFLATE_BUILD_TESTS=OFF"], check=True, env=env)
    subprocess.run(["cmake", "--build", str(build / "codec"), "--parallel", "2", "--target", "libdeflate_static"], check=True, env=env)
    # Seqair's native sources have no public C header. Match main.goose's tag.
    home = Path(env.get("GOOSE_HOME", Path.home() / ".goose"))
    seqair = home / "pkg/github.com/sstadick/goose-seqair/v0.2.1"
    flags = ["-O3", "-march=native"]
    if args.check: flags += ["-g", "-fsanitize=address,undefined", "-fno-sanitize-recover=all", "-fno-omit-frame-pointer"]
    command = ["clang", *flags, str(source), *[str(seqair / "src" / name) for name in ("bases.c", "cram.c", "sam.c")],
               str(package("bgzf.h") / "src/seek.c"), str(package("goosio.h") / "src/goosio.c"),
               str(next((build / "codec").rglob("libdeflate.a"))),
               "-lm", "-pthread", "-ldl", "-o", str(build / ("goose-perbase-check" if args.check else "goose-perbase"))]
    print(" ".join(command), flush=True); subprocess.run(command, check=True, env=env)


if __name__ == "__main__": main()
