struct StringVec <: AbstractVector{String}
    data::Vector{UInt8}
    offsets::Vector{UInt32}

    StringVec() = new(UInt8[], [0])
end

Base.size(x::StringVec) = (length(x.offsets) - 1,)

function Base.getindex(x::StringVec, i::Integer)
    @boundscheck checkbounds(x, i)
    start = @inbounds x.offsets[i] + 1
    stop = @inbounds x.offsets[i + 1]
    return StringView(ImmutableMemoryView(x.data)[start:stop])
end

Base.push!(x::StringVec, str::AbstractString) = push_codeunits!(x, codeunits(str))

function push_codeunits!(x::StringVec, codeunits::AbstractVector{UInt8})
    append!(x.data, codeunits)
    push!(x.offsets, length(x.data))
    return x
end
