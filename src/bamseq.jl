struct BAMSeq <: AbstractString
    data::ImmutableMemoryView{UInt8}
    len::Int
end

Base.ncodeunits(x::BAMSeq) = x.len

const SEQ_LUT = let
    x = zero(UInt128)
    for (i, char) in enumerate("=ACMGRSVTWYHKDBN")
        encoding = UInt8(char) % UInt128
        encoding <<= (8 * i - 8)
        x |= encoding
    end
    x
end

function Base.codeunit(x::BAMSeq, i::Integer)
    i = Int(i)::Int
    @boundscheck checkbounds(1:x.len, i)
    byte = @inbounds x.data[(i + 1) >> 1]
    byte >>= ((4 * isodd(i)) & 7)
    shift = 8 * (byte & 0x0f)
    return (SEQ_LUT >>> shift) % UInt8
end

function Base.thisind(x::BAMSeq, i::Int)
    iszero(i) && return 0
    i == x.len + 1 && return i
    @boundscheck checkbounds(1:x.len, i)
    return i + 1
end

Base.thisind(x::BAMSeq, i::Integer) = thisind(x, Int(i)::Int)

function Base.nextind(x::BAMSeq, i::Int, n::Int = 1)
    @boundscheck checkbounds(0:x.len, i)
    return i + n
end

Base.nextind(x::BAMSeq, i::Integer, n::Integer = 1) = nextind(x, Int(i)::Int, Int(n)::Int)

function Base.iterate(x::BAMSeq, i::Int = 1)
    (i - 1) % UInt ≥ (x.len % UInt) && return nothing
    cu = @inbounds codeunit(x, i)
    return (reinterpret(Char, UInt32(cu) << 24), i + 1)
end

function Base.prevind(x::BAMSeq, i::Int, n::Int = 1)
    @boundscheck checkbounds(1:(x.len + 1), i)
    return i + n
end

Base.prevind(x::BAMSeq, i::Integer, n::Integer = 1) = prevind(x, Int(i)::Int, Int(n)::Int)

function Base.isvalid(x::BAMSeq, i::Integer)
    i = Int(i)::Int
    return in(i, 1:x.len)
end
