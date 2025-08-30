module BioSequencesExt

using BioSequences: LongDNA, BioSequence, NucleicAcidAlphabet, DNA_Gap

using XAM2: XAM2, BAMRecord, has_eq_seq_symbol, HasEq

function XAM2.try_sequence(
        ::Type{LongDNA{4}},
        record::BAMRecord
    )::Union{LongDNA{4}, Nothing, HasEq}
    iszero(record.l_seq) && return nothing
    has_eq_seq_symbol(record) && return HasEq()
    n_seq_bytes = (record.l_seq + 1) >> 1
    offset = record.l_read_name + 4 * record.n_cigar_op
    chunks = Vector{UInt64}(undef, Int(div(record.l_seq + 0x07, 0x08)))
    data = record.data
    GC.@preserve chunks data begin
        src = pointer(data, offset + 1)
        dst = Ptr{UInt8}(pointer(chunks))
        unsafe_copyto!(dst, src, n_seq_bytes)
    end
    seq = LongDNA{4}(chunks, UInt(record.l_seq))
    for i in eachindex(chunks)
        chunk = chunks[i]
        new_chunk = (chunk & 0x0f0f0f0f0f0f0f0f) << 4
        new_chunk |= (chunk & 0xf0f0f0f0f0f0f0f0) >> 4
        chunks[i] = new_chunk
    end
    return seq
end

end # module BioSequencesExt.jl
