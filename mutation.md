# Mutation design

Design notes for record and header mutation in XAM2. Work in progress.

## Guiding principle
Add a mutator when both of these hold:
1. Doing it by hand would break a structural invariant, or the encoding makes it fiddly
   (4-bit sequence packing, shifting the bytes after a field, `0xff` sentinels).
2. There is only one sensible way to do it.

Leave out operations that require a policy choice from the tool, e.g. "trim this read"
(soft or hard clip? What if the trim cuts into an insertion? Update NM/MD?).
Rare needs are served by a text-based escape hatch.
Add new mutators when Overburden or another real user needs them, not speculatively.

## SAM records
`SAMRecord` gets no setters. `BAMRecord` is the only record type that can be edited.
Editing a SAM record means converting it to BAM, editing, and converting back (or writing),
using the existing `try_copy!` in both directions.

* Generic tools (e.g. a read trimmer) work on both SAM and BAM.
* Users never have to build SAM text by hand, which is easy to get wrong
  (phred+33 offset, `*`/`=` conventions, aux number formatting).
* The conversion cost doesn't matter much, because text parsing already dominates
  SAM performance. The small losses from conversion (base case, aux integer widths)
  are acceptable.

`SAMWriter` writes a `SAMRecord`'s stored bytes as they are, and converts a `BAMRecord` with `try_copy!`.

## Headers
Keep the approach of `header_edit.jl`: curated, functional operations returning a new
`SAMHeader`, plus a text escape hatch (edit the text, then reparse with `try_parse_header`).
Headers are small and rarely changed, so performance doesn't matter and immutable rebuilds are fine.

Existing: `add_program`, `add_read_group`, `set_sort_order`, `with_references`.

Worth adding:
* Filter lines with a predicate (e.g. drop `@CO` or certain `@RG` lines)
* Add `@CO` lines

Left out:
* An API for editing individual header tags. Large surface, little gain over editing text.
* Merging headers. This is all policy, e.g. what to do when `@RG` or `@PG` IDs collide.

## BAM records

### Principle: setters never leave a broken record
Every setter call leaves the record consistent and writable.
Constraints between fields are kept either by setting the fields together, or by
making the setter that would break them throw. An order dependency between setters
shows up as an exception, never as corruption. Implicit state, such as the `CG` tag
holding a long CIGAR, is managed by the setters and never exposed.

A setter checks the constraints that involve the fields it changes, against the
resulting state: it never writes a broken state, even if the record was already
broken. Records read from files that are inconsistent can still be edited in
unrelated ways (name, mapq, aux).

### Fixed-length fields
Flags, mapq, placement (refID + pos), next placement (next refID + next pos), tlen.
Setters take `nothing` for "missing", mirroring the getters.
Positions are 1-based, like `left_position` and `next_position`.
Setters throw on bad input, such as a position or reference ID out of range.
References are given either as a 1-based index (like `reference_index`) or a
reference name, since the record holds `refs`.

#### Placement
Reference and position are one concept, *placement*, and are set together.
A record is placed if it has both, and unplaced if it has neither. The valid states:

| State              | `UNMAP` | reference + position             |
|--------------------|---------|----------------------------------|
| Mapped             | unset   | both present                     |
| Unmapped, placed   | set     | both present (usually the mate's)|
| Unmapped, unplaced | set     | both `nothing`                   |

* `set_placement!(record, ref, pos)` and `set_placement!(record, nothing)`.
  `pos` must be within the reference. The alignment's end (pos + the CIGAR's
  reference length) is not checked, since that would tie `set_placement!` and
  `set_cigar!` together, and samtools accepts alignments past the reference end.
  Unplacing throws if `UNMAP` is unset, or if the CIGAR field is a `kSmN` placeholder
  for a long CIGAR (which is only recognised on placed records).
* `set_flags!` throws if the new flags lack `UNMAP` and the record is unplaced.
  So to map an unplaced read, set the placement first, then clear `UNMAP`.
* Mate fields: `set_next_placement!` with the same "both or neither" rule.
  They are not tied to `MUNMAP` or `PAIRED`; those rules are recommended practice
  in the spec, not structural, and belong to fixmate.

### Variable-length fields
Name, CIGAR, seq+qual. Each setter shifts the rest of `data`.
Seq and qual are set together, since they share `l_seq`. `qual = nothing` fills with `0xff`.
Aux data is already editable.

Every setter keeps this true: the record has no CIGAR, or no sequence (`l_seq == 0`),
or the CIGAR's query length equals `l_seq`.
* `set_cigar!(record, cigar)` throws unless `l_seq == 0` or the lengths match.
* `set_sequence!(record, seq, qual)` throws if the record has a CIGAR
  (`n_cigar_op > 0`) and the new length differs from `l_seq`. This only reads a fixed
  field, so the setter never parses the CIGAR. It also rules out changing between a
  missing and a present sequence on a record with a CIGAR.
* A joint setter for CIGAR + seq + qual covers changes in length, e.g. trimming.
  The three fields are adjacent, so it moves the trailing bytes only once.

Since the lengths always match, the long-CIGAR placeholder `kSmN` can always use
`k = l_seq`, as htslib writes it.

Setters never need to parse the existing variable-length content of the record,
except the CIGAR setter's `CG` handling, which uses `aux_data`, not `aux_data_unchecked`.
So they work on records with malformed content, and only throw on bad user input.

### Aliasing
A setter's input may be a view of the record's own data, e.g. setting the name from an aux
value (`set_query_name!(r, aux_data(r)["RX"])`), or `set_cigar!(r, cigar(r))` for a CIGAR
stored in the `CG` tag. Moving the trailing bytes would then corrupt the input.

Setters handle aliasing by views of predictable types (`DataView` in `bam.jl`):
`MemoryView{UInt8}`, and 1-D contiguous byte `SubArray`s of a `Vector`, `Memory` or
`MemoryView`. String inputs are checked by their `codeunits`, which for `StringView`s of
`MemoryView`s and substrings of those is a `MemoryView`. `BAMCIGAR`s are checked by their
`MemoryView`, so CIGAR inputs are always checked. This covers everything the
accessors return, and lets the setters' logic handle only byte arrays.
Other types, e.g. user-defined views, are not checked, so the docstrings of such setters
still warn about aliasing.
`Base.mightalias` is not used, since it is private and unreliable.

The mechanism, without copying:
1. `data_indices(data, x)` compares addresses to find the indices of the bytes of the byte
   array `x` in `data`, or returns `nothing`. For other types, it returns `nothing` at compile time.
2. `replace_region!(write!, record, offset, old_len, new_len, src)` calls `write!(data, src_start)`
   at a point where the input is intact: before moving the trailing bytes when shrinking,
   and after moving them when growing, with `src_start` adjusted if the input was after the
   region. `write!` copies aliased input with a `memmove`, so the input may overlap the region.
3. An input split across the end of the region can't be handled this way when growing.
   Setters must rule it out (for names, the region ends in a null byte, which a valid
   name can't contain), or handle it in an outlined branch, e.g. by copying.

### `bin`
Following the SAM spec, section 5.3: unmapped reads, and reads whose CIGAR consumes no
reference bases, are binned as having length one. Unmapped reads with no position
(-1 in BAM) get `reg2bin(-1, 0) = 4680`.
So `bin` depends on pos, CIGAR and the `UNMAP` flag.
Note: The current `try_copy!(::BAMRecord, ::SAMRecord)` ignores the `UNMAP` flag,
which is a bug.

### Validation and `bin` at write time
Setters keep edited records consistent, but records read from files may not be,
so the writer still checks every record:
`try_write` checks that the CIGAR matches the sequence length, and computes `bin` from
pos + CIGAR + flags. Both need the same pass over the CIGAR ops.

Nothing in XAM2 reads the stored `bin` except the writer (via the fixed-field `memcpy`).
There are two options:
1. *Stateless*: Validate and compute `bin` for every record written.
   The stored `bin` field is then never trusted, or could be dropped.
2. *Stateful*: Only validate records that were edited, tracked by a `geometry_dirty` bit (see below).

The stateful approach is only used if benchmarks show it significantly improves write speed.
The guess is that validating every record is small next to BGZF compression,
and read-only tools like `flagstat` never write.

**Decision: Stateless.** Benchmark on `ecoli.sorted.bam` (4.0M records, 1.1 GB uncompressed),
writing to memory from records already in memory, 8 Julia threads:

| Setting                   | No validation | Validation |
|---------------------------|---------------|------------|
| Default (4 workers, lvl 6)| 2.02 s        | 2.17 s     |
| `compress_level = 1`      | 0.97 s        | 0.97 s     |
| `n_workers = 1`           | 6.77 s        | 6.90 s     |

The validation pass alone (`cigar` + `compute_bin`) takes 0.20 s, about 50 ns per record.
The writer does not check the placement rules, so files that break them still round-trip,
as with samtools.
`compute_bin` reproduces the bins written by samtools for all records in the file.
The stored `bin` field is kept only for the memcpy layout. It is never read,
not set by `try_copy!`, and not part of `==` or `hash`.
A consequence is that records with a malformed CIGAR, or one that doesn't match the
sequence length, can't be written, even if they were read unchanged from a file.

#### The `geometry_dirty` bit (not implemented, since the stateless approach was chosen)
* It is set only by setters that affect `bin` or the CIGAR-vs-sequence check:
  position, CIGAR, seq/qual, and flags when `UNMAP` changes.
  Name, mapq, other flags and aux data do not set it.
* It is placed after the fixed-length fields, so `copy_fixed_fields!` and the memcpy
  layout are unaffected. `copy` must carry it over.
* Initial state:
  * Records from a reader start clean. Malformed input records then pass through
    unchanged (garbage in, garbage out, as samtools does). Reading into the reader's
    record must clear the bit, since with `copy=false` the user may have edited it.
  * `BAMRecord(refs)` starts clean, since it is valid by construction.
  * `try_copy!` from SAM starts clean, since it computes `bin` and checks consistency.
* Writing does *not* clear the bit. Clearing it would turn writing into a mutation of the
  record, which is surprising and a data race if the same record is written from two tasks.
  Recomputing again on a second write is cheap.
* This is a cache invariant, not a structural one, so nothing `@inbounds` depends on it:
  "if not dirty, `bin` is correct for the current pos/CIGAR/flags, or came from the input file".
* The `try_write` docstring must say that only edited records are validated.
  The library only takes responsibility for records it helped create.

### In-place operations
Worth having, because they are fiddly and performance-sensitive:
* A mutable view of the quality scores (for quality binning or masking). Trivial in BAM.
* Setting individual bases in the sequence (`setindex!` on the 4-bit encoding).
* Keeping only a subrange of seq+qual, in place. Shifting nibbles by an odd offset
  is easy to get wrong.

### Left out
* Full read trimming that also adjusts the CIGAR and position. CIGAR clipping helpers,
  if wanted, belong in CIGARStrings.jl.
* MD/NM recomputation, which needs the reference sequence.
* Fixmate and TLEN calculation, and syncing mate fields.
* Records that alias the reader's buffer (read-only records). These would preclude
  mutation, and copy-on-write is not wanted.

## Implementation plan
Each stage is one reviewable commit of at most a few hundred lines.
Stage 2 must come before stages 3-7. Stages 1 and 8 are independent of the rest.
Stage 4b reworks stage 3 code, and must come before stage 5.
The setter naming convention (e.g. `set_mapping_quality!`) is decided in stage 3,
and used by all later stages.

### Stage 1: SAM writer (done)
* `SAMWriter(io, header)`, with `try_write` and `write` for `SAMRecord`.
  A `SAMRecord`'s `data` is the line without the newline, so it's written verbatim.
* `try_write(::SAMWriter, ::BAMRecord)` converts via a scratch `SAMRecord` held by the
  writer, using the existing `try_copy!`.
* Tests: SAM -> SAM and BAM -> SAM round trips on the test files.

### Stage 2: Bin computation, write-time validation, benchmark (done, stateless)
* Extract a `compute_bin(pos, ref_len, flags)` helper from `try_copy_kind!`, following
  the spec for unmapped reads. This fixes the `UNMAP` bug in `try_copy!`.
* Benchmark writing a realistic BAM file, with and without validating and computing
  `bin` for every record (`cigar(record)` + `compute_bin`).
* If the cost is not significant: Implement the stateless approach. In `try_write_kind`,
  return the error from `cigar(record)`, and patch the computed `bin` into the output
  buffer, without writing it back to the record. Consider dropping the stored `bin` field.
* Otherwise, implement the `geometry_dirty` bit as described above:
  * Replace the `@assert sizeof(BAMRecord) == ...` with a `fieldoffset` check,
    since the new field adds padding.
  * Initialize it in `new_uninitialized_record` (`new(data, refs)` leaves it undefined),
    clear it in `reset_fixed_fields!`, `reset!` and `try_read_record!`, carry it in `copy`.
  * Validate and patch `bin` in `try_write_kind` only if the bit is set.
  * Tests set the bit and the position via internal fields, since there are no setters yet.
* Update the `try_write` docstring accordingly.

### Stage 3: Fixed-field setters (done)
Flags, mapq, pos, next pos, reference, next reference, tlen.
`nothing` for missing, 1-based positions and reference indices, reference names accepted.
Throw on out-of-range input, leaving the record unchanged.

Conventions, used by all later stages:
* Setters are named `set_<getter>!`, e.g. `set_mapping_quality!`, and return the record.
  The reference setters are `set_reference!` and `set_next_reference!`.
* A setter accepts exactly the values its getter can return, so
  `getter(setter!(record, x)) == x`. E.g. mapq 255 and tlen 0 throw; use `nothing`.
* Positions are not checked against the reference length or the CIGAR, so setters
  can be called in any order.

### Stage 4: Shifting helper and name setter (done)
* Internal `replace_region!(record, offset, old_len, new_len)`: check the block-size
  invariant before modifying anything, then resize `data` and move the trailing bytes.
* `set_query_name!` validates the name against the spec's QNAME regex, which also ensures at
  most 254 bytes and no null bytes (invariant 3). `nothing` is stored as `*`, and `"*"`
  itself throws, following the stage 3 convention.
* Aliasing by the record's own name or aux values is handled, see "Aliasing".
* Names are handled as bytes only for string types known to be UTF-8 (`UTF8String`:
  `String`, `StringView`, and substrings of those). Other strings may have wider code
  units (UTF-16), or ASCII-range code units that are not ASCII characters (UTF-7), so
  `to_utf8_string` converts them with `String` in an outlined method, and is a no-op otherwise.

### Stage 4b: Placement (reworks stage 3) (done)
* Replace `set_reference!` + `set_left_position!` with `set_placement!`, and
  `set_next_reference!` + `set_next_position!` with `set_next_placement!`.
* Unplacing throws if `UNMAP` is unset, or the CIGAR field is a placeholder.
* `set_flags!` throws if the new flags lack `UNMAP` and the record is unplaced.

### Stage 5: CIGAR setter (done)
* Accept `Union{Nothing, BAMCIGAR}`. Users convert from `CIGAR` themselves, so they
  opt into the cost of encoding.
* Throw if the CIGAR is empty (the getter returns `nothing` for no CIGAR), or if
  `l_seq != 0` and the query length differs from `l_seq`.
* More than 65535 ops: store the placeholder `kSmN` (`k = l_seq`, `m` = reference
  length) and the CIGAR in the `CG` tag, like `try_copy!`. Throw if the record is
  unplaced (the placeholder would not be recognised), if `k` or `m` don't fit in 28 bits,
  or if the aux data is malformed (a `CG` tag in it would not be read).
* Always delete an existing `CG` tag, since it is stale after setting a new CIGAR.
  If the aux data is malformed and the new CIGAR is short, skip this: `load_aux_cigar`
  ignores malformed aux data, so the tag is inert.
* Aliasing: replace the CIGAR field with `replace_region!`, then delete or rewrite
  `CG`. Input that would be split by the move, and aliased input for a long CIGAR
  (e.g. `set_cigar!(r, cigar(r))` with the CIGAR in `CG`), is copied in an outlined branch.
  Since `BAMCIGAR`s are always checked, the docstring needs no aliasing warning.

### Stage 6: Sequence and quality setters (done)
* `set_sequence!(record, seq, qual)` and `set_cigar_and_sequence!(record, cigar, seq, qual)`.
* Seq is a string (or `nothing`) of `=ACMGRSVTWYHKDBN`, the characters BAM can store.
  Lowercase is accepted and stored as uppercase. Unlike when parsing SAM (where, like
  htslib, other letters and `.` become `N`), other characters throw, since in a
  setter they are more likely a bug. A `BAMSeq` (e.g. from
  another record) is copied without re-encoding. Empty strings throw, like empty CIGARs.
* Qual is PHRED scores without offset, like `quality_scores`; `nothing` fills with
  `0xff`. Lengths must match, and a missing sequence requires missing quality.
  Scores must be in `0:93`, so the record can be written as SAM.
* `set_sequence!` throws if the record has a CIGAR and the length changes.
  The joint setter moves the trailing bytes once, and handles long CIGARs.
* Input is validated before shifting `data`. Aliased input (e.g. trimming with the
  record's own seq and qual) is copied in an outlined branch; efficient in-place
  trimming is stage 7b. Qual of types other than `DataView`s is not checked for
  aliasing, so the docstrings warn about it.
* Later, if needed: a fast path for `LongDNA{4}` in the BioSequences extension. Its
  4-bit encoding matches BAM's for the nucleotide codes.

### Stage 7a: Simple in-place operations
A mutable view of the quality scores, and setting individual bases.
These change no lengths.

### Stage 7b: Keeping a subrange of seq+qual
In place, handling odd nibble offsets.
Test exhaustively against a naive implementation (decode, slice, re-encode),
for all subranges of all sequences up to ~9 bases.

### Stage 8: Header additions
Filter lines with a predicate, and add `@CO` lines, in the style of `header_edit.jl`.

### Afterwards
* With the setters, `BAMRecord(refs)` plus setters covers building records from scratch
  (`TODO.md` items 3 and 4). A keyword constructor can be added later if Overburden needs it.
* Update `TODO.md`: Item 3's builder idea is superseded by this design.
