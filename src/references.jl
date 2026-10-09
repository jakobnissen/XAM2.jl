"""
    Reference

A reference sequence, as stored in a [`References`](@ref).
Has the public properties `name::AbstractString` and `length::Int`.

See also: [`References`](@ref)
"""
struct Reference
    name::StringView{ImmutableMemoryView{UInt8}}
    length::Int
end

Base.:(==)(a::Reference, b::Reference) = (a.name == b.name) & (a.length == b.length)
Base.hash(x::Reference, h::UInt) = hash(x.length, hash(x.name, hash(Reference, h)))

function Base.show(io::IO, x::Reference)
    print(io, "Reference(")
    show(io, x.name)
    return print(io, ", ", x.length, ')')
end

"""
    References <: AbstractVector{Reference}

Immutable list of the reference sequences of a BAM file, e.g. chromosomes.
Every reference has a unique, non-empty name and a length in `1:typemax(Int32)`.

Readers and writers share their `References` with every record they create,
such that records can look up their reference names and lengths.
Use [`references`](@ref) to obtain the `References` of a reader or writer.

Construct with `References(names, lengths)`, or with `References(header::SAMHeader)`
to use the `@SQ` lines of a header. These constructors throw an `ArgumentError`
if the references are invalid.

See also: [`Reference`](@ref)

# Examples
```jldoctest
julia> refs = References(["chr1", "chr2"], [1000, 500]);

julia> refs[2]
Reference("chr2", 500)

julia> refs[1].length
1000
```
"""
const ReferenceIndex = Dict{StringView{ImmutableMemoryView{UInt8}}, Int32}

# Holds a lazily built `ReferenceIndex`. It's atomic because `References`
# may be shared between threads.
mutable struct LazyReferenceIndex
    @atomic index::Union{Nothing, ReferenceIndex}
end

struct References <: AbstractVector{Reference}
    names::StringVec
    lengths::Vector{Int32}
    # Name => zero-based reference ID, as stored in BAM records.
    # This is built on the first lookup by name, because it takes a significant amount
    # of memory for references with many sequences, and is not needed for reading.
    lazy_index::LazyReferenceIndex

    global function unsafe_new_references(names::StringVec, lengths::Vector{Int32})
        return new(names, lengths, LazyReferenceIndex(nothing))
    end
end

"""
    references(x::Union{BAMReader, BAMWriter})::References

Get the reference sequences of `x`.
"""
function references end

Base.size(x::References) = (length(x.lengths),)
Base.IndexStyle(::Type{References}) = IndexLinear()

function Base.getindex(x::References, i::Integer)
    @boundscheck checkbounds(x, i)
    return Reference(@inbounds(x.names[i]), @inbounds(x.lengths[i]) % Int)
end

# Return the zero-based reference ID of `name`, or `nothing` if it is not present.
function reference_id(refs::References, name::AbstractString)::Union{Nothing, Int32}
    return get(get_index(refs), name, nothing)
end

function get_index(refs::References)::ReferenceIndex
    lazy = refs.lazy_index
    index = @atomic :acquire lazy.index
    index === nothing || return index
    # If multiple threads race to build the index, they each build an identical
    # index, so it doesn't matter which one is stored.
    index = ReferenceIndex()
    sizehint!(index, length(refs))
    for i in eachindex(refs.names)
        index[@inbounds(refs.names[i])] = (i - 1) % Int32
    end
    @atomic :release lazy.index = index
    return index
end

# Validate the references. Return a `References`, or an `XAMErrorKind` on
# invalid references.
function try_build_references(
        names::StringVec,
        lengths::Vector{Int32},
    )::Union{References, XAMErrorKind}
    length(names) == length(lengths) || return XAMErrors.invalid_reference
    length(lengths) > typemax(Int32) && return XAMErrors.invalid_reference
    # Only used to check for duplicates, so it's freed after construction.
    seen = Set{StringView{ImmutableMemoryView{UInt8}}}()
    sizehint!(seen, length(names))
    for i in eachindex(names, lengths)
        name = @inbounds names[i]
        # The name is stored null-terminated in BAM, so it can't contain null bytes.
        # Names longer than this can't be stored in BAM, since the stored length
        # includes the terminating null byte.
        if isempty(name) || ncodeunits(name) ≥ typemax(Int32) || in(0x00, codeunits(name))
            return XAMErrors.invalid_reference
        end
        @inbounds(lengths[i]) < 1 && return XAMErrors.invalid_reference
        n_before = length(seen)
        push!(seen, name)
        length(seen) == n_before && return XAMErrors.duplicate_reference
    end
    return unsafe_new_references(names, lengths)
end

function References(
        names::AbstractVector{<:AbstractString},
        lengths::AbstractVector{<:Integer},
    )
    if length(names) != length(lengths)
        throw(ArgumentError("Number of reference names and lengths must be the same"))
    end
    stringvec = StringVec()
    for name in names
        push!(stringvec, name)
    end
    reflens = Vector{Int32}(undef, length(lengths))
    for (i, len) in zip(eachindex(reflens), lengths)
        if len < 1 || len > typemax(Int32)
            throw(ArgumentError("Reference length must be in 1:2147483647"))
        end
        reflens[i] = len % Int32
    end
    return unwrap_references(try_build_references(stringvec, reflens))
end

function References(header::SAMHeader)
    stringvec = StringVec()
    reflens = Vector{Int32}(undef, length(header.seqs))
    for (i, sq) in enumerate(header.seqs)
        push_codeunits!(stringvec, sq.sn)
        if iszero(sq.ln) || sq.ln > typemax(Int32) % UInt32
            throw(ArgumentError("Reference length must be in 1:2147483647"))
        end
        reflens[i] = sq.ln % Int32
    end
    return unwrap_references(try_build_references(stringvec, reflens))
end

function unwrap_references(x::Union{References, XAMErrorKind})::References
    x isa References && return x
    msg = if x == XAMErrors.duplicate_reference
        "Reference names must be unique"
    else
        "Reference names must be non-empty, shorter than 2147483647 bytes, " *
            "and not contain null bytes, and there must be fewer than 2147483648 references"
    end
    throw(ArgumentError(msg))
end

# Check that the @SQ lines of the header, if any, agree with `refs`.
function check_header_references(header::SAMHeader, refs::References)::Nothing
    seqs = header.seqs
    isempty(seqs) && return nothing
    if length(seqs) != length(refs)
        throw(ArgumentError("Number of @SQ lines in header does not match number of references"))
    end
    for (i, sq) in enumerate(seqs)
        ref = @inbounds refs[i]
        if sq.sn != codeunits(ref.name) || sq.ln != ref.length
            throw(ArgumentError("@SQ line $(i) in header does not match reference $(i)"))
        end
    end
    return nothing
end

# Sentinel values in a reference translation map. Non-negative values
# are reference IDs of the destination.
const REFMAP_MISSING = Int32(-2)
const REFMAP_LENGTH_MISMATCH = Int32(-3)

# Build the map from reference IDs in `src` to IDs in `dst`, where the map
# at index `i` contains the destination ID of source ID `i - 1`.
# Return `nothing` if `src` and `dst` have identical content, such that IDs
# need no translation.
function build_reference_map(src::References, dst::References)::Union{Nothing, Vector{Int32}}
    # Check content equality first, which avoids building the index of `dst`
    is_identical = src.lengths == dst.lengths &&
        src.names.offsets == dst.names.offsets &&
        src.names.data == dst.names.data
    is_identical && return nothing
    dst_index = get_index(dst)
    map = Vector{Int32}(undef, length(src))
    for i in eachindex(map, src.lengths)
        dst_id = get(dst_index, @inbounds(src.names[i]), nothing)
        @inbounds map[i] = if dst_id === nothing
            REFMAP_MISSING
        elseif @inbounds(dst.lengths[dst_id + 1]) != @inbounds(src.lengths[i])
            # Same name with a different length is almost certainly a different
            # assembly, so we refuse to translate it.
            REFMAP_LENGTH_MISMATCH
        else
            dst_id
        end
    end
    return map
end

# Translate the reference ID `id` (where -1 means missing) using a map
# from `build_reference_map`.
function translate_reference_id(map::Vector{Int32}, id::Int32)::Union{Int32, XAMErrorKind}
    id == -1 && return id
    (id < -1 || id ≥ length(map)) && return XAMErrors.invalid_reference_id
    new_id = @inbounds map[id + 1]
    new_id == REFMAP_MISSING && return XAMErrors.unknown_reference
    new_id == REFMAP_LENGTH_MISMATCH && return XAMErrors.reference_length_mismatch
    return new_id
end
