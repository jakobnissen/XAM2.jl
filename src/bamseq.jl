"""
    BAMSeqString <: AbstractString

Lazy, immutable view of the 4-bit encoded sequence of a [`BAMRecord`](@ref), returned by
[`sequence(::BAMRecord)`](@ref sequence). Bases are decoded to ASCII on access.

Its code units are `UInt8`, and it contains only the ASCII characters
`"=ACMGRSVTWYHKDBN"`, so its characters, code units and indices correspond one to one.

Passing a `BAMSeqString` to [`set_sequence!`](@ref) copies its encoded bases without
re-encoding them.

!!! warning
    A `BAMSeqString` aliases its record. Mutating the record will invalidate it.
"""
struct BAMSeqString <: AbstractString
    # Invariant: `length(data) == cld(len, 2)`
    data::ImmutableMemoryView{UInt8}
    # Each byte stores two characters, so we need this field
    # to distinguish an odd from an even length
    len::Int
end

Base.ncodeunits(x::BAMSeqString) = x.len

const SEQ_LUT = Tuple(b"=ACMGRSVTWYHKDBN")

# Like htslib, letters that are not IUPAC codes, and `.`, are encoded as `N`.
# Other bytes are invalid.
const ENCODE_LUT = let
    result = fill(0xff, 128)
    for char in ['A':'Z'; 'a':'z'; '.']
        result[Int(char) + 1] = 0x0f
    end
    for (i, char) in enumerate("=ACMGRSVTWYHKDBN")
        encoding = i - 1
        for c in [char, lowercase(char)]
            result[Int(c) + 1] = encoding
        end
    end
    Tuple(result)
end

# Like `ENCODE_LUT`, but only the 16 BAM codes (and their lowercase) are valid.
# Used by setters, where other letters are more likely a bug than an intended `N`.
const STRICT_ENCODE_LUT = let
    result = fill(0xff, 128)
    for (i, char) in enumerate("=ACMGRSVTWYHKDBN")
        encoding = i - 1
        for c in [char, lowercase(char)]
            result[Int(c) + 1] = encoding
        end
    end
    Tuple(result)
end

Base.codeunit(::BAMSeqString) = UInt8

function Base.codeunit(x::BAMSeqString, i::Integer)
    i = Int(i)::Int
    @boundscheck checkbounds(1:x.len, i)
    byte = @inbounds x.data[(i + 1) >> 1]
    byte >>= ((4 * isodd(i)) & 7)
    return SEQ_LUT[(byte & 0x0f) + 1]
end

function Base.thisind(x::BAMSeqString, i::Int)
    iszero(i) && return 0
    i == x.len + 1 && return i
    @boundscheck checkbounds(1:x.len, i)
    return i
end

Base.thisind(x::BAMSeqString, i::Integer) = thisind(x, Int(i)::Int)

function Base.nextind(x::BAMSeqString, i::Int, n::Int = 1)
    @boundscheck checkbounds(0:x.len, i)
    return i + n
end

Base.nextind(x::BAMSeqString, i::Integer, n::Integer = 1) = nextind(x, Int(i)::Int, Int(n)::Int)

function Base.iterate(x::BAMSeqString, i::Int = 1)
    (i - 1) % UInt ≥ (x.len % UInt) && return nothing
    cu = @inbounds codeunit(x, i)
    return (reinterpret(Char, UInt32(cu) << 24), i + 1)
end

function Base.prevind(x::BAMSeqString, i::Int, n::Int = 1)
    @boundscheck checkbounds(1:(x.len + 1), i)
    return i - n
end

Base.prevind(x::BAMSeqString, i::Integer, n::Integer = 1) = prevind(x, Int(i)::Int, Int(n)::Int)

function Base.isvalid(x::BAMSeqString, i::Integer)
    i = Int(i)::Int
    return in(i, 1:x.len)
end

function bamseq!(dst::MutableMemoryView{UInt8}, src::ImmutableMemoryView{UInt8})
    target_len = (length(src) + 1) >>> 1
    if length(dst) < target_len
        throw(DimensionMismatch("dst size must be cld(length(src), 2)"))
    end
    for i in 1:target_len
        @inbounds dst[i] = 0x00
    end
    for (i, byte) in enumerate(src)
        shift = isodd(i) ? 4 : 0
        byte > 127 && @goto bad_encoding
        encoding = @inbounds ENCODE_LUT[byte + 1]
        encoding == 0xff && @goto bad_encoding
        @inbounds dst[(i + 1) >>> 1] |= encoding << (shift & 7)
    end
    return BAMSeqString(ImmutableMemoryView(dst[1:target_len]), length(src))
    @label bad_encoding
    throw(ArgumentError("BAM sequence can only contain letters, `=` and `.`"))
end
