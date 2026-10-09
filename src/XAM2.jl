module XAM2

using BufferIO: BufferIO,
    AbstractBufReader,
    AbstractBufWriter,
    BufReader,
    VecWriter,
    get_nonempty_buffer,
    get_unflushed,
    consume,
    fill_buffer,
    get_buffer,
    grow_buffer,
    line_views

using MemoryViews:
    MemoryViews,
    ImmutableMemoryView,
    MemoryView,
    MutableMemoryView,
    split_last,
    split_each,
    split_at

using StringViews: StringView
using BGZFLib: BGZFReader, BGZFWriter, virtual_position, VirtualOffset
using CIGARStrings: CIGARStrings, AbstractCIGAR, BAMCIGAR, CIGAR, CIGARError, unsafe_switch_memory
using XAMAuxData: XAMAuxData, AuxTag, try_auxtag

using Automa: Automa, @re_str, onenter!, onexit!, onall!, rep, rep1, opt

export BAMRecord,
    BAMReader,
    BAMWriter,
    SAMRecord,
    SAMReader,
    SAMHeader,
    SAMFlag,
    SAMFlags,
    Flag,
    XAMRecord,
    push,
    StringVec,
    References,
    references,
    get_header,
    query_name,
    flags,
    reference_name,
    reference_length,
    next_reference_name,
    left_position,
    sequence,
    quality_string,
    quality_scores,
    aux_data,
    is_mapped,
    is_first,
    is_last,
    mapping_span,
    XAMError,
    XAMErrorKind,
    XAMErrors

public try_parse!,
    try_write,
    Reference,
    try_parse_header,
    HasEq,
    SAMHeaderError

"""
    module XAMErrors

This module is used as a namespace for the [`XAMErrorKind`](@ref) enum.
"""
module XAMErrors

    """
        XAMErrorKind

    Enum that contains information about a [`XAMError`](@ref).
    New variants may be added in future, minor releases.
    The integer value of these kinds are not stable over minor releases.

    julia> kind = XAMErrors.invalid_cigar;

    julia> throw(XAMError(kind))
    ERROR:
    [...]
    """
    @enum XAMErrorKind::UInt8 begin
        truncated_file
        invalid_block_length
        closed_stream
        header_truncated
        header_bad_magic_bytes
        header_top_bit_set
        invalid_cigar
        missing_fields # not enough tab fields
        invalid_integer
        overflow
        line_too_long # more than 4 GiB line
        qname_too_long
        qual_no_seq
        qual_seq_mismatch
        invalid_reference # reference with invalid name or length
        duplicate_reference # two references with the same name
        invalid_reference_id # record has a reference ID outside the references
        unknown_reference # record's reference is not present in the destination
        reference_length_mismatch # same reference name, but different lengths
    end
end # module XAMErrors

using .XAMErrors: XAMErrorKind

"""
    XAMError <: Exception

An exception thrown when reading or writing invalid SAM or BAM data.

Instances of this error has a public `.kind` property of type
`XAMErrorKind`. You can use this programmatically check the kind
of exception.
"""
struct XAMError <: Exception
    kind::XAMErrorKind
end

@noinline throw_xam_error(kind::XAMErrorKind) = throw(XAMError(kind))

function Base.showerror(io::IO, error::XAMError)
    kind = error.kind
    s = if kind == XAMErrors.truncated_file
        "Input data too is truncated; does not contain whole record or header"
    elseif kind == XAMErrors.invalid_block_length
        unreachable()
    end
    return print(io, s)
end

struct Uninit end

"""
    HasEq()

Zero-sized struct indicating that a BAM sequence contains `=`, and therefore cannot
be returned as a `LongSequence`.
"""
struct HasEq end

const CTX = Automa.CodeGenContext(; generator = :goto, clean = true)

@static if isdefined(Base, :memoryindex)
    memory_index(x::MemoryRef) = Base.memoryindex(x)
else
    memory_index(x::MemoryRef) = Core.memoryrefoffset(x)
end

# Index in underlying memory of first index of MemoryView or Vector
memory_index(v::Vector) = memory_index(Base.cconvert(Ptr, v))
memory_index(m::MemoryView) = first(only(parentindices(m)))

function tryparse_base10(::Type{T}, mem::ImmutableMemoryView{UInt8}) where {T <: Union{UInt8, UInt16, UInt32, UInt64}}
    isempty(mem) && return nothing
    n = UInt64(0)
    has_overflowed = false
    bad_byte = false
    for i in mem
        bad_byte |= i ∉ 0x30:0x39
        (n, f) = Base.Checked.mul_with_overflow(n, UInt64(0x0a))
        has_overflowed |= f
        (n, f) = Base.Checked.add_with_overflow(n, (i - 0x30) % UInt64)
        has_overflowed |= f
    end
    ((n > typemax(T)) | has_overflowed | bad_byte) && return nothing
    return n % T
end

function tryparse_base10(::Type{T}, mem::ImmutableMemoryView{UInt8}) where {T <: Union{Int8, Int16, Int32}}
    isempty(mem) && return nothing
    start = @inbounds mem[1]
    (is_positive, maxval, mem) = if start == UInt8('-')
        length(mem) < 2 && return nothing
        (false, abs(Int64(typemin(T))) % UInt64, @inbounds(mem[2:end]))
    else
        (true, typemax(T) % UInt64, mem)
    end
    n = @inline tryparse_base10(UInt64, mem)
    isnothing(n) && return nothing
    n > maxval && return nothing
    n = n % T
    return is_positive ? n : -n
end

@inline function unsafe_load_integer_le(::Type{T}, mem::ImmutableMemoryView{UInt8}, i::Int) where {T <: Integer}
    x = zero(T)
    shift = 0
    mask = 8 * sizeof(T) - 1
    for offset in 0:(sizeof(T) - 1)
        byte = @inbounds mem[i + offset]
        x |= (byte % T) << (shift & mask)
        shift += 8
    end
    return x
end

function inbounds_write_base10!(mem::MutableMemoryView{UInt8}, n::Unsigned)::Int
    T = typeof(n)
    if iszero(n)
        @inbounds mem[1] = 0x30
        return 1
    end
    n_digits = 0
    @inbounds while !iszero(n)
        n_digits += 1
        (n, rem) = divrem(n, T(10))
        mem[n_digits] = (rem % UInt8) + 0x30
    end
    rv = n_digits
    @inbounds for i in 1:(n_digits >> 1)
        (mem[i], mem[rv]) = (mem[rv], mem[i])
        rv -= 1
    end
    return n_digits
end

function reg2bin(span::UnitRange{Int32})::UInt16
    # TODO: remove these checks, also the one at the bottom
    (start, stop) = (first(span), last(span))
    @assert start in Int32(0):Int32(2^29 - 1)
    @assert stop in (start):(Int32(2^29 - 1))

    # If start is zero, the read is unmapped, and so the range
    # must be 0:0 exactly
    @assert !iszero(start) || iszero(stop)

    # Convert start to an offset to make it zero-indexed.
    # Note that this will be -1 when span is 0:0.
    # Then this will use arithmetic shift right.
    start -= one(start)
    stop -= one(stop)

    result = if start >> 14 == stop >> 14
        div((1 << 15) - 1, 7) + (start >> 14)
    elseif start >> 17 == stop >> 17
        div((1 << 12) - 1, 7) + (start >> 17)
    elseif start >> 20 == stop >> 20
        div((1 << 9) - 1, 7) + (start >> 20)
    elseif start >> 23 == stop >> 23
        div((1 << 6) - 1, 7) + (start >> 23)
    elseif start >> 26 == stop >> 26
        div((1 << 3) - 1, 7) + (start >> 26)
    else
        0
    end

    @assert result in 0:Int(typemax(UInt16))

    return result % UInt16
end

function reg2bins!(bins::Vector{UInt16}, span::UnitRange{Int32})
    (start, stop) = (first(span) - Int32(1), last(span) - Int32(1))
    # @assert start in Int32(1):Int32(2^29 - 1)
    # @assert stop in start:(Int32(2^29 - 1))

    # If start is zero, the read is unmapped, and so the range
    # must be 0:0 exactly
    # @assert !iszero(start) || iszero(stop)

    # TODO: Check the 0,0 case
    idxs = (1, 9, 73, 585, 4681)
    shifts = (26, 23, 20, 17, 14)

    ranges = map(1:5) do i
        idx = @inbounds(idxs[i]) % UInt16
        shift = @inbounds shifts[i] & 63
        (idx + (start >> shift) % UInt16):(idx + (stop >> shift) % UInt16)
    end

    resize!(bins, sum(length, ranges; init = 0) + 1)

    @inbounds bins[1] = 0

    i = 1
    for range in ranges, k in range
        i += 1
        @inbounds bins[i] = k % UInt16
    end

    return bins
end

# This internal function _assumes_ that:
# * If old_cigar is a AbstractCIGAR, it aliases `old_vector`
# * old_vector == new_vector
function copy_cigar!(
        old_cigar::Union{Nothing, Uninit, CIGARError, AbstractCIGAR},
        old_vector::Vector{UInt8},
        new_vector::Vector{UInt8},
    )
    if !isa(old_cigar, AbstractCIGAR)
        return old_cigar
    end
    # First: Which index in `old_vector` does the backing memory of
    # `old_cigar` start?
    cigar_mem = MemoryView(old_cigar)
    index_cigar_mem = memory_index(cigar_mem)
    index_old_vec = memory_index(old_vector)
    cigar_mem_start = index_cigar_mem - index_old_vec + 1

    # Now, copy the bytes of cigar_mem into `new_vector` as those same
    # indices
    # TODO: Inbounds once we've tested correctness
    dst_mem = MemoryView(new_vector)[cigar_mem_start:((cigar_mem_start + length(cigar_mem)) - 1)]
    copy!(dst_mem, cigar_mem) # TODO: Inbounds once we've checked
    dst_mem = ImmutableMemoryView(dst_mem)
    new_cigar = unsafe_switch_memory(old_cigar, dst_mem)
    return new_cigar
end

"""
    seek_first_record(reader::Union{SAMReader, BAMReader}) -> reader

Seek `reader` back to the first record.
"""
function seek_first_record end

include("flags.jl")
include("stringvec.jl")
include("header.jl")
include("references.jl")
include("bamseq.jl")
include("sam.jl")
include("bam.jl")

"""
Alias for Union{SAMRecord, BAMRecord}
"""
const XAMRecord = Union{SAMRecord, BAMRecord}

"""
Alias for Union{SAMReader, BAMReader}
"""
const XAMReader = Union{SAMReader, BAMReader}

"""
    query_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the query name (QNAME) of `record`, or `nothing` if the data
is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

!!! warning
    The returned string aliases the record. Mutating the record will
    invalidate the string.
"""
function query_name end

"""
    reference_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the reference name (RNAME) of `record`, or `nothing` if the data
is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

For `SAMRecord`, this function does _not_ validate that the reference
name is present in the SAM header.

!!! warning
    For `SAMRecord`, this aliases the record, whereas for `BAMRecord`,
    the data aliases the `References` of the record.
    Mutating the underlying record will invalidate the string if returned
    from a `SAMRecord`
"""
function reference_name end

"""
    next_reference_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the reference name of the next segment in the template (RNEXT) of `record`,
or `nothing` if the data is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

!!! warning
    For `SAMRecord`, this aliases the record, whereas for `BAMRecord`,
    the data aliases the `References` of the record.
    Mutating the underlying record will invalidate the string if returned
    from a `SAMRecord`
"""
function next_reference_name end

"""
    reference_length(record::BAMRecord)::Union{Nothing, Int}

Return the length of the reference that `record` is mapped to, or `nothing`
if the record is not mapped.
"""
function reference_length end

"""
    left_position(record::XAMAuxData)::Union{Nothing, Int}

Return the 1-based leftmost (i.e. 5' of the reference) mapped position
of `record`, or `nothing` if the data is not available.
"""
function left_position end

"""
    cigar(record::SAMRecord)::Union{Nothing, CIGAR, CIGARError}
    cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError}

Get the `cigar` of `record` as a `CIGARStrings.CIGAR` or
`CIGARStrings.BAMCIGAR`, or return `nothing` if the data is not available,
or a `CIGARStrings.CIGARError` if the cigar is malformed.

This operation will validate the underlying data, and then cache
the constructed cigar such that subsequent calls to this function
is almost free.

!!! warning
    The returned `AbstractCIGAR` object aliases the record.
    Mutating the record will invalidate the cigar.
"""
function cigar end

"""
    flags(record::XAMRecord)::SAMFlags

Get the flags of `record`.
"""
flags(x::XAMRecord) = x.flags

"""
    is_mapped(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.UNMAP`.
"""
is_mapped(x::XAMRecord) = !in(UNMAP, flags(x))

"""
    is_first(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.READ1`.
"""
is_first(x::XAMRecord) = in(READ1, flags(x))

"""
    is_last(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.READ2`.
"""
is_last(x::XAMRecord) = in(READ2, flags(x))

"""
    get_header(::XAMReader)::SAMHeader

Get the header of the reader.
"""
function get_header end

"""
    mapping_span(record::XAMRecord)::Union{Nothing, UnitRange{Int}}

Return the 1-based positions in the reference that `record` is mapped to,
or `nothing` if this data is not available.

E.g. a record with CIGAR `10M100S` mapping to reference position 5
would map to `5:14`.

This data requires the record to contain a valid `POS` and `CIGAR` fields.
"""
function mapping_span(record::XAMRecord)
    left = @something left_position(record) return nothing
    cig = cigar(record)
    if cig === nothing || cig isa CIGARError
        return nothing
    end
    return left:(left + CIGARStrings.ref_length(cig) - 1)
end

"""
    try_sequence(::Type{LongDNA{4}}, rec::BAMRecord)::Union{Nothing, HasEq, LongDNA{4}}

Get the sequence of `rec` as a `LongDNA{4}` biosequence.

This function returns `nothing` if the sequence is missing (i.e. if `sequence(rec)`
would return `nothing`)

This function returns `XAM2.HasEq()` if the sequence contains a `=`. In this case,
use the method `sequence(::LongDNA{4}, ::BAMRecord, ::BioSequence{<:DNAAlphabet})`.

!!! note
    This function requires `BioSequences` to be loaded
"""
function try_sequence end

end # module XAM2
