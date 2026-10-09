module BioSequencesExt

using BioSequences: LongDNA, BioSequence, DNAAlphabet, DNA_Gap, isgap

using XAM2: XAM2, BAMRecord, first_eq_seq_pos, HasEq, cigar
using CIGARStrings: pos2pos, query, ref, CIGARError

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
        rec::BAMRecord,
        ref::BioSequence{<:DNAAlphabet}
    )::Union{Nothing, LongDNA{4}}

Return the read sequence of `rec` as a `LongDNA{4}` biosequence.

This function is complicated by the fact that the sequence can be stored as `=`, which indicates
the same base as the equivalent reference position. Therefore, if `=` is present in the seq, this
function relies on mapping information to get the corresponding reference base.
    
This function returns `nothing` in any of the below cases:
* The sequence information of the record is not present
* `=` is present, but the record does not have a reference
* `=` is present, and either the POS or CIGAR is missing from the record

This function throws:
* If the passed `ref` length does not match `rec`'s reference length, throw an `ArgumentError`
* If `=` is present, and the record has an invalid CIGAR, throw an `ErrorException`
"""
function XAM2.sequence(
        ::Type{LongDNA{4}},
        record::BAMRecord,
        ref::BioSequence{<:DNAAlphabet}
    )
    iszero(record.l_seq) && return nothing # no seq
    seq = get_seq(record)
    eq_pos = first_eq_seq_pos(record)
    eq_pos === nothing && return seq
    reflen = XAM2.reference_length(record)
    reflen === nothing && return nothing # no ref
    length(ref) == reflen || throw(ArgumentError("Given reference has wrong length"))
    span = XAM2.mapping_span(record)
    span === nothing && return nothing
    cig = cigar(record)
    cig === nothing && return nothing # no CIGAR
    cig isa Union{CIGARStrings.CIGARError, XAM2.XAMErrorKind} && error("Record has invalid CIGAR")
    refview = view(ref, span)
    gap_positions = (i for i in eq_pos:length(seq) if isgap(seq[i]))
    for (gap_position, translation) in pos_to_pos(query, ref, cig, gap_positions)
        @assert translation.kind === CIGARStrings.pos
        seq[gap_position] = refview[translation.pos]
    end
    return seq
end

function get_seq(record::BAMRecord)
    n_seq_bytes = (record.l_seq + 1) >> 1
    offset = record.l_read_name + 4 * record.n_cigar_op
    chunks = Vector{UInt64}(undef, Int(div(record.l_seq + 0x07, 0x08)))
    data = record.data
    # Copy over data
    GC.@preserve chunks data begin
        src = pointer(data, offset + 1)
        dst = Ptr{UInt8}(pointer(chunks))
        unsafe_copyto!(dst, src, n_seq_bytes)
    end
    # Reverse each nuc in representation. This relies on
    # internal repr of LongSequence, but that is unlikely to change.
    for i in eachindex(chunks)
        chunk = chunks[i]
        new_chunk = (chunk & 0x0f0f0f0f0f0f0f0f) << 4
        new_chunk |= (chunk & 0xf0f0f0f0f0f0f0f0) >> 4
        chunks[i] = new_chunk
    end
    return LongDNA{4}(chunks, UInt(record.l_seq))
end

end # module BioSequencesExt.jl
