#!/usr/bin/env bash
set -euo pipefail
exec python3 "$(dirname "${GOOSE_DOC:-../goose-doc/main.goose}")/doctest.py" \
  --cc clang --symbols api.txt --preamble tests/doc_preamble.goose --require-examples \
  "$@" perbase.goose
