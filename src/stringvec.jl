struct StringVec <: AbstractVector{String}
    data::Vector{UInt8}
    offsets::Vector{UInt32}

    StringVec() = new(UInt8[], [0])
end

Base.size(x::StringVec) = (length(x.offsets) - 1,)

function Base.getindex(x::StringVec, i::Integer)
    @boundscheck checkbounds(x, i)
    start = x.offsets[i] + 1
    stop = x.offsets[i + 1]
    return StringView(ImmutableMemoryView(x.data)[start:stop])
end

function Base.push!(x::StringVec, str::AbstractString)
    cu = codeunits(str)::AbstractVector{UInt8}
    append!(x.data, cu)
    push!(x.offsets, length(x.data))
    return x
end