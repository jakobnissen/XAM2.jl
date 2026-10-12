"""
    QualityView <: AbstractVector{UInt8}

Lazy view of the quality scores of a record, which converts between the two
encodings of quality scores on access:
SAM stores PHRED scores plus 33, as printable ASCII characters, and BAM stores the PHRED
scores themselves.

A `QualityView` is returned by [`quality_scores(::SAMRecord)`](@ref quality_scores), where
it subtracts 33 from each byte of the SAM record, and its bytes back the string returned by
[`quality_string(::BAMRecord)`](@ref quality_string), where it adds 33 to each byte of the
BAM record. The arithmetic wraps around, so out-of-range bytes never throw.

The elements are computed on every access. When the scores are accessed more than once,
or when a contiguous array is needed (e.g. for `MemoryView` or `pointer`),
materialize the view with `collect`, which returns a `Vector{UInt8}`. `collect`, and
`copyto!` to other arrays, are efficient.

!!! warning
    A `QualityView` aliases its record. Mutating the record will invalidate it.
    In particular, the destination of `copyto!` must not share memory with the view,
    e.g. be the record's data, since `copyto!` does not check for aliasing.

# Examples
```jldoctest
julia> record = XAM2.parse!(SAMRecord(), MemoryView("r\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACG\\t!+I"));

julia> scores = quality_scores(record)
3-element XAM2.QualityView:
 0x00
 0x0a
 0x28

julia> collect(scores)
3-element Vector{UInt8}:
 0x00
 0x0a
 0x28
```
"""
struct QualityView <: AbstractVector{UInt8}
    data::ImmutableMemoryView{UInt8}
    # Added to each byte of `data`, wrapping around
    offset::UInt8
end

# Converts quality scores from SAM to BAM encoding
sam_to_bam_quality(data::ImmutableMemoryView{UInt8}) = QualityView(data, -UInt8(33))

# Converts quality scores from BAM to SAM encoding
bam_to_sam_quality(data::ImmutableMemoryView{UInt8}) = QualityView(data, UInt8(33))

Base.size(x::QualityView) = (length(x.data),)
Base.IndexStyle(::Type{QualityView}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(x::QualityView, i::Int)
    @boundscheck checkbounds(x, i)
    return @inbounds(x.data[i]) + x.offset
end

# The generic `copyto!` checks for aliasing with `Base.mightalias`, which is slow for this
# type, since it falls back to hashing with `objectid`. It also can't detect aliasing, since
# `MemoryView` doesn't support it. So these methods don't check for aliasing, and their
# destination must not alias the view, as documented.
Base.copyto!(dst::AbstractArray, src::QualityView) = copyto!(dst, firstindex(dst), src, 1, length(src))

# Resolves a method ambiguity with Base
Base.copyto!(dst::PermutedDimsArray, src::QualityView) = copyto!(dst, firstindex(dst), src, 1, length(src))

function Base.copyto!(dst::AbstractArray, dstart::Integer, src::QualityView, sstart::Integer, n::Integer)
    dstart = Int(dstart)::Int
    sstart = Int(sstart)::Int
    n = Int(n)::Int
    n < 0 && throw(ArgumentError("Number of elements to copy must be non-negative"))
    iszero(n) && return dst
    checkbounds(dst, dstart:(dstart + n - 1))
    checkbounds(src, sstart:(sstart + n - 1))
    @inbounds for i in 0:(n - 1)
        dst[dstart + i] = src[sstart + i]
    end
    return dst
end
