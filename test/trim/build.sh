#!/usr/bin/env bash
# Trim smoke test: Compile app.jl with `juliac --trim=safe`, and run it.
# Requires Julia >= 1.12. Usage: test/trim/build.sh
# Set JULIAC_ENV to reuse an environment with JuliaC installed.
set -euo pipefail
cd "$(dirname "$0")"

JULIAC_ENV="${JULIAC_ENV:-$(mktemp -d)}"
julia --project="$JULIAC_ENV" -e 'using Pkg; Pkg.add("JuliaC")'
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# --bundle fails if the target directory already exists.
rm -rf build
# Note: If verification fails, juliac prints "Verifier error" lines, and exits non-zero.
julia --project="$JULIAC_ENV" -e 'using JuliaC; JuliaC._main_cli(ARGS)' -- \
    --output-exe xam2trim --project=. --trim=safe --bundle build app.jl

OUT="$(mktemp -d)/out.bam"
./build/bin/xam2trim ../data/example.bam "$OUT" ../data/example.sam
