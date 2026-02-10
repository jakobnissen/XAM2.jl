# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

XAM2.jl is a high-performance Julia library for reading and writing SAM/BAM bioinformatics alignment files. It is a from-scratch reimplementation focused on speed, using low-level memory manipulation, Automa.jl-based parsing, and BGZF compression.

This project is under development, so some code reflects work-in-progress and temporary workarounds as opposed to the result of active design descisions. 

## Development Commands

```bash
# Run tests (no test files exist yet)
JULIA_TEST_FAILFAST=true julia --startup=no --project -e 'using Pkg; Pkg.test()'

# Load the package in a REPL
julia --startup=no --project -e 'using XAM2'

# Instantiate dependencies (many are local path dependencies)
julia --startup=no --project -e 'using Pkg; Pkg.instantiate()'
```

## Key Dependencies (local path)

Several dependencies are developed alongside this package and referenced by local path in `Project.toml` (`[sources]` section):
- `Automa.jl` - Finite state machine code generation for parsing
- `BGZFLib.jl` - BGZF (blocked gzip) compression for BAM files
- `CIGARStrings.jl` - CIGAR string parsing and manipulation
- `MemoryViews.jl` - Zero-copy memory views
- `XAMAuxData.jl` - SAM/BAM auxiliary data handling

These are expected at `../Automa.jl`, `../BGZFLib.jl`, etc. relative to this repo.

## Architecture

### Record Types
- **`BAMRecord`** (`bam.jl`): Mutable struct holding BAM alignment data. Fixed-length fields are stored as struct fields (loaded via `memcpy` from raw bytes). Variable-length data (qname, cigar, seq, qual, aux) stored in a `Vector{UInt8}`. Each record holds a back-reference to its reader/writer for resolving reference IDs to names.
- **`SAMRecord`** (`sam.jl`): Mutable struct holding SAM alignment data. The raw line is stored in `data::Vector{UInt8}`, with offset fields pointing into it. Parsed using Automa-generated state machine code.

### Reader/Writer Types
- **`BAMReader`/`BAMWriter`** (`bam.jl`): Wrap `BGZFReader`/`BGZFWriter`. The reader reuses a single internal `BAMRecord` across iterations (controlled by `copy` kwarg). Header and reference sequence metadata are parsed eagerly on construction.
- **`SAMReader`** (`sam.jl`): Wraps an `AbstractBufReader`. Same copy-on-iterate pattern as BAMReader.

### Parsing Strategy
- SAM parsing uses Automa.jl to generate a `goto`-based state machine (`SAM_MACHINE` with `SAM_ACTIONS` dict). The `parse_record!` function is `@eval`-generated from the machine.
- BAM parsing reads fixed-length fields via `memcpy` into struct fields, then reads variable-length data into the record's vector.
- SAM/BAM headers are parsed by a separate Automa machine (`HEADER_MACHINE` in `header.jl`).

### Performance Patterns
- Fixed struct fields are loaded/stored via `memcpy` + `pointer_from_objref` at computed offsets (`FIXED_BAM_OFFSET`, `FIXED_SAM_OFFSET`), bypassing Julia's normal field access.
- Records are reused during iteration to avoid allocation; `copy` creates an independent snapshot.
- `Uninit` sentinel type is used for lazy CIGAR parsing (computed on first access, then cached).
- Integer parsing (`tryparse_base10`) is hand-rolled for performance with overflow checking.

### Supporting Types
- **`SAMFlags`** (`flags.jl`): Bitfield `AbstractSet{Flag}` wrapping `UInt16`, with set operations.
- **`BAMSeq`** (`bamseq.jl`): `AbstractString` that decodes 4-bit-per-base BAM sequence encoding on the fly.
- **`StringVec`** (`stringvec.jl`): Compact vector of strings backed by a single `Vector{UInt8}` with offset table.
- **`SAMHeader`** (`header.jl`): Parsed header with HD, SQ, RG, PG, and CO lines.

### Extension
- `ext/BioSequencesExt.jl`: Package extension loaded when `BioSequences` is available. Provides `try_sequence` and `sequence` methods to convert BAM sequences to `LongDNA{4}`, handling the `=` (reference-equal) encoding.

## Error Handling

Errors are represented as `XAMError` enum values (not exceptions) returned from parsing functions like `try_read!`. The iterator (`Base.iterate`) converts these to thrown errors. Reader objects store the first error in an `error` field and refuse further reads after an error.

## Julia Version

Requires Julia >= 1.11.6.

## When writing commits
Do not add yourself (the agent) as co-author
