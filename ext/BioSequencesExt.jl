module BioSequencesExt

using BioSequences: BioSequences, LongDNA, BioSequence, DNAAlphabet, isgap

using XAM2: XAM2,
    BAMRecord,
    HasEq,
    XAMError,
    XAMErrors,
    cigar,
    left_position,
    reference_length,
    first_eq_seq_pos,
    seq_offset,
    unsafe_sequence_bytes

using MemoryViews: MemoryView, MutableMemoryView

using CIGARStrings: CIGARStrings, pos_to_pos, query, ref, ref_length

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
    )::Union{Nothing, LongDNA{4}, XAMError}

Return the read sequence of `record` as a `LongDNA{4}` biosequence, where `reference`
is the sequence of the record's reference.

A BAM sequence may contain `=`, which means the base is the same as the reference
base it is aligned to. Such bases are looked up in `reference`, using the position
and CIGAR of the record.

Return `nothing` if:
* The record has no sequence
* The sequence contains `=`, and the record has no reference, position or CIGAR

Return an error value if the sequence contains `=`, and:
* The CIGAR is malformed: Return the `XAMError` from [`XAM2.cigar`](@ref)
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
    )::Union{Nothing, LongDNA{4}, XAMError}
    iszero(record.l_seq) && return nothing
    seq = get_seq(record)
    eq_pos = @something first_eq_seq_pos(record) return seq
    reflen = @something reference_length(record) return nothing
    length(reference) == reflen || throw(ArgumentError("Given reference has wrong length"))
    left = @something left_position(record) return nothing
    cig = @something cigar(record) return nothing
    cig isa XAMError && return cig
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
    return get_seq(unsafe_sequence_bytes(record), record.l_seq % Int)
end

# Decode the packed BAM sequence `bytes` of length `l_seq` to a `LongDNA{4}`.
# `bytes` must have `cld(l_seq, 2)` bytes.
function get_seq(bytes::MemoryView{UInt8}, l_seq::Int)::LongDNA{4}
    n_seq_bytes = cld(l_seq, 2)
    length(bytes) == n_seq_bytes || throw(ArgumentError("Wrong number of sequence bytes"))
    # A LongDNA{4} stores 16 bases per chunk, the first base in the lowest bits.
    # The unused bits of the last chunk must be zero.
    chunks = zeros(UInt64, cld(l_seq, 16))
    GC.@preserve chunks bytes begin
        unsafe_copyto!(Ptr{UInt8}(pointer(chunks)), pointer(bytes), n_seq_bytes)
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

# A mutable view of the sequence of a `BAMRecord`.
# The 4-bit encodings of BAM and `DNAAlphabet{4}` are the same, except that the code 0
# is `=` in BAM, and a gap in BioSequences.
# Invariant: `length(data) == cld(len, 2)`, and `len > 0`.
# The sequence has no `=` when the view is created, and the view never writes one,
# so it has no gaps, unless the record is modified by other means.
struct BAMBioSequence <: BioSequence{DNAAlphabet{4}}
    data::MutableMemoryView{UInt8}
    len::Int
end

"""
    try_sequence(::Type{BioSequence}, record::BAMRecord)::Union{Nothing, HasEq, BioSequence{DNAAlphabet{4}}}

Return a mutable view of the sequence of `record` as a `BioSequence{DNAAlphabet{4}}`.
Like [`sequence(record)`](@ref sequence), this does not copy, unlike the methods that
return a `LongDNA{4}`. Setting an element of the view sets the base in the record,
e.g. `seq[5] = DNA_N`. The length of the view can't change.

This function returns `nothing` if the sequence is missing (i.e. if `sequence(record)`
would return `nothing`).

This function returns `XAM2.HasEq()` if the sequence contains a `=`, which means the
base is the same as the reference base, and which `DNAAlphabet{4}` can't represent.
In this case, use the method
`sequence(::Type{LongDNA{4}}, ::BAMRecord, ::BioSequence{<:DNAAlphabet})`.
Checking for `=` takes time proportional to the length of the sequence.

Since a `=` can't be represented, the view never contains gaps. Setting an element to
a gap throws an `ArgumentError`. To store a `=`, use [`XAM2.unsafe_sequence_bytes`](@ref).

The type of the view is an implementation detail. It can't be sliced: To get a
subsequence, copy it to a `LongDNA{4}` first with `copy` or `LongDNA{4}`, which is fast.

!!! warning
    The returned view aliases the record. Any change to the record that moves its
    sequence, e.g. setting its name or CIGAR, invalidates the view, after which writing
    to it corrupts the record.

!!! note
    This function requires `BioSequences` to be loaded

See also: [`XAM2.unsafe_sequence_bytes`](@ref), [`quality_scores`](@ref)
"""
function XAM2.try_sequence(::Type{BioSequence}, record::BAMRecord)::Union{Nothing, HasEq, BAMBioSequence}
    bytes = unsafe_sequence_bytes(record)
    isempty(bytes) && return nothing
    first_eq_seq_pos(record) === nothing || return HasEq()
    return BAMBioSequence(bytes, record.l_seq % Int)
end

Base.length(x::BAMBioSequence) = x.len
BioSequences.encoded_data_eltype(::Type{BAMBioSequence}) = UInt

# Like for `LongSequence`, these are called after bounds checking, e.g. by `getindex`.
# BAM stores the first base of each byte in the high nibble.
@inline function BioSequences.extract_encoded_element(x::BAMBioSequence, i::Integer)
    i = Int(i)::Int
    byte = @inbounds x.data[(i + 1) >>> 1]
    return UInt((byte >>> ((4 * isodd(i)) & 7)) & 0x0f)
end

@inline function BioSequences.encoded_setindex!(x::BAMBioSequence, v::UInt, i::Integer)
    # The code of a gap is 0, which is `=` in BAM
    iszero(v) && throw_gap_error()
    i = Int(i)::Int
    shift = (4 * isodd(i)) & 7
    index = (i + 1) >>> 1
    old = @inbounds x.data[index]
    @inbounds x.data[index] = (old & ~(0x0f << shift)) | ((v % UInt8) << shift)
    return x
end

@noinline function throw_gap_error()
    throw(ArgumentError("Cannot set a gap in a BAM sequence, since it would be stored as `=`"))
end

# The view can't be created or resized, so new sequences are `LongDNA{4}`s.
# Converting copies the packed bytes, instead of the generic element-wise conversion.
LongDNA{4}(x::BAMBioSequence) = get_seq(x.data, x.len)
Base.copy(x::BAMBioSequence) = LongDNA{4}(x)
Base.similar(::BAMBioSequence, len::Integer) = LongDNA{4}(undef, len)

end # module BioSequencesExt
