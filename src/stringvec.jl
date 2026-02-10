"""
    StringVec <: AbstractVector{StringView{ImmutableMemoryView{UInt8}}}

Compact vector of strings, backed by contiguous memory.
Beside the core `AbstractVector` interface, this type supports:

* `Base.push!(::StringVec, ::AbstractString)`

```jldoctest
julia> sv = StringVec()
0-element StringVec

julia> push!(sv, "Hello!")
1-element StringVec:
 "Hello!"
```
"""
struct StringVec <: AbstractVector{StringView{ImmutableMemoryView{UInt8}}}
    data::Vector{UInt8}

    # Has length(StringVec) + 1 elements. E.g. when content is "a", "", "abc",
    # this is [0, 1, 1, 4]
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
