module XAM2

using BufferIO: BufferIO,
    AbstractBufReader,
    AbstractBufWriter,
    BufReader,
    CursorReader,
    VecWriter,
    get_nonempty_buffer,
    get_unflushed,
    consume,
    fill_buffer,
    get_buffer,
    line_views

using MemoryViews:
    ImmutableMemoryView,
    MemoryView,
    MutableMemoryView,
    split_last,
    split_each,
    split_at

using StringViews: StringView
using BGZFLib: BGZFReader, BGZFWriter, virtual_position, virtual_seek, VirtualOffset
using CIGARStrings: CIGARStrings, AbstractCIGAR, BAMCIGAR, CIGAR, CIGARError, unsafe_switch_memory
using XAMAuxData: XAMAuxData, AuxTag, try_auxtag

using Automa: Automa, @re_str, onenter!, onexit!, onall!, rep, rep1, opt

export BAMRecord,
    BAMReader,
    BAMWriter,
    SAMRecord,
    SAMReader,
    SAMHeader

module XAMErrors

    @enum XAMError::UInt8 begin
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

    end
end # module XAMErrors

using .XAMErrors: XAMError

struct Uninit end

struct HasEq end

const CTX = Automa.CodeGenContext(; generator = :goto, clean = true)

function to_memory(v)
    result = Memory{eltype(v)}(undef, length(v))
    copyto!(result, v)
    return result
end


@static if isdefined(Base, :memoryindex)
    memory_index(x::MemoryRef) = memoryindex(x)
else
    memory_index(x::MemoryRef) = Core.memoryrefoffset(x)
end

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

    ranges = map((1, 9, 73, 585, 4681)) do idx
        idx = idx % UInt16
        (idx + (start >> 26) % UInt16):(idx + (stop >> 26) % UInt16)
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


copy_cigar(x::Union{Nothing, Uninit, CIGARError}, _::Vector{UInt8}, _::Vector{UInt8}) = x

function copy_cigar(old_cigar::AbstractCIGAR, old_v::Vector{UInt8}, cp::Vector{UInt8})
    # Find the view of `cp` that corresponds to the same data as the
    # view that backs `old_cigar`
    cigarmem = MemoryView(old_cigar)
    cigar_offset = memory_index(cigarmem)
    v_offset = memory_index(old_v)
    cigar_v_offset = cigar_offset - v_offset
    new_mem = ImmutableMemoryView(cp)[(cigar_v_offset + 1):(cigar_v_offset + length(cigarmem))]
    return unsafe_switch_memory(old_cigar, new_mem)
end

"""
    seek_first_record(reader::Union{SAMReader, BAMReader}) -> reader

Seek `reader` back to the first record.
"""
function seek_first_record end


include("flags.jl")
include("stringvec.jl")
include("header.jl")
include("bamseq.jl")
include("sam.jl")
include("bam.jl")

const XAMRecord = Union{SAMRecord, BAMRecord}

flags(x::XAMRecord) = x.flags
is_mapped(x::XAMRecord) = !in(UNMAP, flags(x))
is_forward(x::XAMRecord) = in(READ1, flags(x))
is_reverse(x::XAMRecord) = in(READ2, flags(x))

function mapping_span(record::XAMRecord)
    left = @something left_position(record) return nothing
    cig = cigar(record)
    if cig === nothing || cig isa CIGARError
        return nothing
    end
    return left:(left + CIGARStrings.ref_length(cig) - 1)
end


# TODO: Does this make sense? We already have ref information elsewhere.
# maybe delete it.
# TODO: Add more information?
# TODO: Add a show method
# This should be super cheap to construct - make it cheaper!
struct Reference
    name::StringView{ImmutableMemoryView{UInt8}}
    len::Int
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
