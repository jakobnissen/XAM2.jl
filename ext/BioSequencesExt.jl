module BioSequencesExt

using BioSequences: LongDNA, BioSequence, DNAAlphabet, isgap

using XAM2: XAM2,
    BAMRecord,
    HasEq,
    XAMError,
    XAMErrors,
    cigar,
    left_position,
    reference_length,
    first_eq_seq_pos,
    seq_offset

using CIGARStrings: CIGARStrings, CIGARError, pos_to_pos, query, ref, ref_length

function XAM2.try_sequence(
        ::Type{LongDNA{4}},
        record::BAMRecord
    )::Union{LongDNA{4}, Nothing, HasEq}
    iszero(record.l_seq) && return nothing
    first_eq_seq_pos(record) === nothing || return HasEq()
    return get_seq(record)
end

"""
    sequence(
        ::Type{LongDNA{4}},
        record::BAMRecord,
        reference::BioSequence{<:DNAAlphabet}
    )::Union{Nothing, LongDNA{4}, CIGARError, XAMError}

Return the read sequence of `record` as a `LongDNA{4}` biosequence, where `reference`
is the sequence of the record's reference.

A BAM sequence may contain `=`, which means the base is the same as the reference
base it is aligned to. Such bases are looked up in `reference`, using the position
and CIGAR of the record.

Return `nothing` if:
* The record has no sequence
* The sequence contains `=`, and the record has no reference, position or CIGAR

Return an error value if the sequence contains `=`, and:
* The CIGAR is malformed: Return the `CIGARError` or `XAMError` from [`XAM2.cigar`](@ref)
* The alignment extends past the end of `reference`: Return an `XAMError` of kind
  `XAMErrors.invalid_position`
* A `=` is not aligned to a reference base, e.g. it is soft clipped or inserted:
  Return an `XAMError` of kind `XAMErrors.invalid_sequence`

Throw an `ArgumentError` if the length of `reference` differs from the length of the
record's reference.
"""
function XAM2.sequence(
        ::Type{LongDNA{4}},
        record::BAMRecord,
        reference::BioSequence{<:DNAAlphabet}
    )::Union{Nothing, LongDNA{4}, CIGARError, XAMError}
    iszero(record.l_seq) && return nothing
    seq = get_seq(record)
    eq_pos = @something first_eq_seq_pos(record) return seq
    reflen = @something reference_length(record) return nothing
    length(reference) == reflen || throw(ArgumentError("Given reference has wrong length"))
    left = @something left_position(record) return nothing
    cig = @something cigar(record) return nothing
    cig isa Union{CIGARError, XAMError} && return cig
    # The alignment's end is not validated when reading
    if left + ref_length(cig) - 1 > length(reference)
        return XAMError(XAMErrors.invalid_position)
    end
    # `=` is decoded as a gap, since both are encoded as 0
    gap_positions = (i for i in eq_pos:length(seq) if isgap(seq[i]))
    # The CIGAR's query length is the sequence length, else `cigar` returns an error.
    # Translated reference positions are relative to the start of the alignment.
    for (gap_position, translation) in pos_to_pos(query, ref, cig, gap_positions)
        translation.kind == CIGARStrings.pos || return XAMError(XAMErrors.invalid_sequence)
        seq[gap_position] = reference[left + translation.pos - 1]
    end
    return seq
end

function get_seq(record::BAMRecord)::LongDNA{4}
    l_seq = record.l_seq % Int
    n_seq_bytes = cld(l_seq, 2)
    # A LongDNA{4} stores 16 bases per chunk, the first base in the lowest bits.
    # The unused bits of the last chunk must be zero.
    chunks = zeros(UInt64, cld(l_seq, 16))
    data = record.data
    # Copy over data. Invariant: The sequence fits in `data`.
    GC.@preserve chunks data begin
        src = pointer(data, seq_offset(record) + 1)
        dst = Ptr{UInt8}(pointer(chunks))
        unsafe_copyto!(dst, src, n_seq_bytes)
    end
    # BAM stores the first base of each byte in the high nibble, so swap the nibbles.
    # The 4-bit encodings of BAM and BioSequences are the same.
    for i in eachindex(chunks)
        chunk = chunks[i]
        new_chunk = (chunk & 0x0f0f0f0f0f0f0f0f) << 4
        new_chunk |= (chunk & 0xf0f0f0f0f0f0f0f0) >> 4
        chunks[i] = new_chunk
    end
    # If the length is odd, the last byte has a padding nibble, which may be nonzero
    n_used_bits = (4 * l_seq) % 64
    if !iszero(n_used_bits)
        chunks[end] &= (UInt64(1) << n_used_bits) - 1
    end
    return LongDNA{4}(chunks, UInt(l_seq))
end

end # module BioSequencesExt
