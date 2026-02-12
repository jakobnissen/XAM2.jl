module SAMHeaderErrors
    @enum SAMHeaderErrorKind::UInt8 begin
        invalid_sort_order
        invalid_topology
        no_vn_field
        invalid_vn_field
        invalid_ln_field
        no_sn_field
        no_ln_field
        no_id_field
        duplicate_tag
    end

end

using .SAMHeaderErrors: SAMHeaderErrorKind

"""
    SAMHeaderError <: Exception

Exception thrown, or returned, from parsing SAM headers.
This type is public, but its fields are currently internal.
"""
struct SAMHeaderError <: Exception
    line::Int
    # kind or parser error failure
    kind::Union{SAMHeaderErrorKind, Tuple{Int, UInt8}}
end

function Base.showerror(io::IO, error::SAMHeaderError)
    kind = error.kind
    if kind isa Tuple{Int, UInt8}
        (col, byte) = kind
        return print(
            io,
            "Cannot parse SAM, got unexpected byte ",
            repr(byte),
            " at line ",
            error.line,
            " column ",
            col,
            "."
        )
    else
        s = if kind == SAMHeaderErrors.invalid_sort_order
            "In HD line, SO field is not \"unknown\", \"unsorted\", \"queryname\" or \"coordinate\""
        elseif kind == SAMHeaderErrors.invalid_topology
            "TP field in SQ line must be \"linear\" or \"circular\""
        elseif kind == SAMHeaderErrors.no_vn_field
            "In HD line, did not find mandatory VN field"
        elseif kind == SAMHeaderErrors.invalid_ln_field
            "Could not parse LN field in SQ header as 32-bit integer"
        elseif kind == SAMHeaderErrors.no_ln_field
            "SQ line does not have mandatory LN field"
        elseif kind == SAMHeaderErrors.no_sn_field
            "SQ line does not have mandatory SN field"
        elseif kind == SAMHeaderErrors.no_id_field
            "RG or PG header line does not have mandatory ID field"
        elseif kind == SAMHeaderErrors.duplicate_tag
            "Duplicate tags found in header line"
        else
            unreachable()
        end
        prefix = "In SAM header, "
        if error.line isa Int
            prefix *= "in line $(error.line), "
        end
        return print(io, prefix * s)
    end
end

const HEADER_MACHINE = let
    lf = onenter!(re"\n", :newline)
    newline = opt(re"\r") * lf
    header_field = onexit!(onenter!(re"[ -~]+", :mark), :field)
    tag = onexit!(re"[A-Za-z][A-Za-z0-9]", :tag)
    header_segment = onexit!(tag * re":" * header_field, :segment)
    header_content = header_segment * rep(re"\t" * header_segment)
    onexit!(onenter!(header_content, :mark_content), :header_content)

    header_start = onenter!(re"@", :header)
    hd = onexit!(header_start * re"HD\t" * header_content, [:check_unique, :hd])

    # CO can contain arbitrary UTF8 data, but not a newline. I don't bother validating
    # whether it's valid UTF8
    co = onexit!(header_start * onenter!(re"CO\t[^\n]*", :enter_co), :co)
    sq = onexit!(header_start * re"SQ\t" * header_content, [:check_unique, :sq])
    rg = onexit!(header_start * re"RG\t" * header_content, [:check_unique, :rg])
    pg = onexit!(header_start * re"PG\t" * header_content, [:check_unique, :pg])

    # This is slightly complex, because the newline after a CO field must be \n, not
    # \r?\n. If it was the latter, then observing a \r would make the state machine
    # ambiguous.
    non_hd = co | sq | rg | pg
    non_hd_newline = (co * lf * rep(newline)) | ((sq | rg | pg) * rep1(newline))

    header = (
        re"" |
            (hd * rep(newline)) |
            opt(hd * newline) * (non_hd | (rep(non_hd_newline) * non_hd))
    )
    Automa.compile(header)
end

const HEADER_ACTIONS = Dict{Symbol, Expr}(
    :newline => :(line_number += 1),
    :mark => :(mark = p),
    :mark_content => :(content_mark = p),
    :field => :(field = @inbounds(data[mark:(p - 1)])),
    :tag => quote
        tag = (@inbounds(data[p - 2]) % UInt16) | ((@inbounds(data[p - 1]) % UInt16) << 8)
        push!(tags, tag)
    end,
    :segment => :(push!(fields, HeaderField(tag, field))),
    :header => quote
        empty!(tags)
        empty!(fields)
    end,
    :header_content => :(field_mem = @inbounds(data[content_mark:(p - 1)])),
    :check_unique => quote
        if has_duplicate_tags(tags)
            return SAMHeaderError(line_number, SAMHeaderErrors.duplicate_tag)
        end
    end,
    :hd => quote
        hd = try_parse_hd(fields, field_mem)
        hd isa SAMHeaderError && return hd
    end,
    :enter_co => :(content_mark = p + 3),
    :co => quote
        comment_data = @inbounds(data[content_mark:(p - 1)])

        # Since we allow \r in comments, if the line ends with \r\n, the \r
        # needs to be trimmed away
        if !isempty(comment_data) && @inbounds(comment_data[end]) == UInt8('\r')
            comment_data = @inbounds comment_data[1:(end - 1)]
        end
        push!(cos, comment_data)
    end,
    :sq => quote
        sq = try_parse_sq(line_number, fields, field_mem)
        sq isa SAMHeaderError && return sq
        push!(sqs, sq)
    end,
    :rg => quote
        id = try_parse_rg_pg(line_number, fields)
        id isa SAMHeaderError && return id
        push!(rgs, HeaderRG(field_mem, id, length(fields)))
    end,
    :pg => quote
        id = try_parse_rg_pg(line_number, fields)
        id isa SAMHeaderError && return id
        push!(pgs, HeaderPG(field_mem, id, length(fields)))
    end,
)

@noinline function throw_property_error(x, sym)
    s = "Type $(typeof(x)) has no such property: $(sym)"
    error(s)
end

struct HeaderError <: Exception
    line::Int
    str::String
end

struct HeaderField
    tag::UInt16 # the two-byte tag stored in little-endian format
    data::ImmutableMemoryView{UInt8}
end

function str_as_u16(s)
    ncodeunits(s) == 2 || error()
    return (@inbounds(codeunit(s, 1)) % UInt16) | ((@inbounds(codeunit(s, 2)) % UInt16) << 8)
end

"""
    module SortOrders

Namespace module for `SortOrder`
"""
module SortOrders

"""
    SortOrder

1-byte enum with the available sort orders in the SAM/BAM format.
Currently, as per the SAM specs, the following are defined:
`unknown`, `unsorted`, `queryname`, `corrdinate` 
"""
@enum SortOrder::UInt8 unknown unsorted queryname coordinate
export SortOrder, unknown, unsorted, queryname, coordinate
end

using .SortOrders: SortOrders, SortOrder

function parse_mem(
        ::Type{SortOrder},
        mem::ImmutableMemoryView{UInt8}
    )::Union{SortOrder, SAMHeaderErrorKind}
    return if mem == b"unknown"
        SortOrders.unknown
    elseif mem == b"unsorted"
        SortOrders.unsorted
    elseif mem == b"queryname"
        SortOrders.queryname
    elseif mem == b"coordinate"
        SortOrders.coordinate
    else
        return SAMHeaderErrors.invalid_sort_order
    end
end

"""
    module Topologies

Namespace module for `Topology`
"""
module Topologies

"""
    Topology

1-byte enum with the available topologies in the SAM/BAM format.
Currently, as per the SAM specs, the following are defined:
`linear`, `circular`.
If not specified, `linear` is assumed. 
"""
@enum Topology::UInt8 linear circular
export Topology, linear, circular
end

using .Topologies: Topologies, Topology

"""
    AbstractHeader <: AbstractDict{AuxTag, ImmutableMemoryView{UInt8}}

A SAM header line, e.g. `@SQ\\tSN:sequence1\\tLN:2015`.
This is an immutable `AbstractDict` where the keys are `XAMAuxData.AuxTag`, and
the values are an `ImmutableMemoryView{UInt8}` with the corresponding value.

The content of the line, not including the first `@`-prefixed tag and the first
tab, can be obtained with `MemoryView(::AbstractHeader)`.

See also: [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> sam_header = get_header(SAMReader(IOBuffer(sam_data)));

julia> sq_header = sam_header.seqs[1]
XAM2.HeaderSQ with 2 entries:
  AuxTag("SN") => [0x43, 0x48, 0x52, 0x4f, 0x4d, 0x4f, 0x53, 0x4f, 0x4d, 0x45, 0x5f, 0x49]
  AuxTag("LN") => [0x31, 0x30, 0x30, 0x39, 0x38, 0x30, 0x30]

julia> sq_header["LN"] |> String
"1009800"

julia> MemoryView(sq_header) |> String
"SN:CHROMOSOME_I\\tLN:1009800"
```
"""
abstract type AbstractHeader <: AbstractDict{AuxTag, ImmutableMemoryView{UInt8}} end

function Base.getproperty(x::AbstractHeader, symbol::Symbol)
    if symbol ∈ propertynames(x)
        return getfield(x, symbol)
    end
    throw_property_error(x, symbol)
end

MemoryViews.MemoryView(x::AbstractHeader) = getfield(x, :fields)

"""
    HeaderHD <: AbstractHeader

A struct representing a `@HD` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `major::UInt32`: Major version of the "VN" field
* `minor::UInt32`: Minor version of the "VN" field
* `sort_order::SortOrder`. Defaults to `SortOrders.unknown`

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> hd_header = get_header(SAMReader(IOBuffer(sam_data))).header;

julia> hd_header.major
0x00000001

julia> hd_header.sort_order == SortOrders.unsorted
true

julia> hd_header["VN"] |> String
"1.6"
```
"""
struct HeaderHD <: AbstractHeader
    # All content of the header are stored in contiguous memory
    fields::ImmutableMemoryView{UInt8}
    len::Int # number of fields
    major::UInt32
    minor::UInt32
    sort_order::SortOrder
end

Base.propertynames(::HeaderHD) = (:major, :minor, :sort_order)

"""
    HeaderSQ <: AbstractHeader

A struct representing a `@SQ` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `sn::ImmutableMemoryView{UInt8}`: Name of sequence
* `ln::UInt32`: Length of the sequence
* `topology::Topology`. Defaults to `Topologies.linear`

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> sq_header = get_header(SAMReader(IOBuffer(sam_data))).seqs[1];

julia> sq_header.topology == Topologies.linear
true

julia> sq_header.ln
0x000f6888

julia> sq_header.sn |> String
"CHROMOSOME_I"
```
"""
struct HeaderSQ <: AbstractHeader
    # All content of the header are stored in contiguous memory
    fields::ImmutableMemoryView{UInt8}
    sn::ImmutableMemoryView{UInt8} # a view into `fields`
    len::Int # number of fields
    ln::UInt32
    topology::Topology # linear by default
end

Base.propertynames(::HeaderSQ) = (:sn, :ln, :topology)

"""
    HeaderRG <: AbstractHeader

A struct representing a `@RG` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `id::ImmutableMemoryView{UInt8}`: Read group identifier

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)
"""
struct HeaderRG <: AbstractHeader
    # All content of the header are stored in contiguous memory
    fields::ImmutableMemoryView{UInt8}
    id::ImmutableMemoryView{UInt8} # a view into `field`
    len::Int # number of fields
end

Base.propertynames(::HeaderRG) = (:id,)

"""
    HeaderPG <: AbstractHeader

A struct representing a `@PG` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `id::ImmutableMemoryView{UInt8}`: Read group identifier

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> pg_header = get_header(SAMReader(IOBuffer(sam_data))).programs[1];

julia> pg_header.id |> String
"samtools"

julia> MemoryView(pg_header) |> String
"ID:samtools\\tPN:samtools\\tVN:1.23\\tCL:samtools view -h ex1.bam"
```
"""
struct HeaderPG <: AbstractHeader
    # Same fields as HeaderRG
    fields::ImmutableMemoryView{UInt8}
    id::ImmutableMemoryView{UInt8}
    len::Int
end

Base.propertynames(::HeaderPG) = (:id,)

"""
    SAMHeader

Object representing a SAM or BAM header, i.e. all the lines in a SAM
file starting with `@`.

Obtain a `SAMHeader` from a [`SAMReader`])(ref) or [`BAMReader`](@ref)
using [`get_header`](@ref).

Obtain the header in SAM format with `MemoryView(::SAMHeader)`.


A `SAMHeader` has the public properties listed below. More may be added in
future releases. Within each `Vector` of `AbstractHeader`, the elements
appear in the same order as in the source SAM/BAM file. 
* `header::Union{Nothing, HeaderHD}`
* `seqs::Vector{HeaderSQ}`
* `comments::Vector{ImmutableMemoryView{UInt8}}`: This represents all the `@CO`
  lines, without the leading `"@CO\\t"`.
* `read_groups::Vector{HeaderRG}`
* `programs::Vector{HeaderPG}`

See also: [`AbstractHeader`](@ref), [`HeaderHD`](@ref)

# Examples
```jldoctest
julia> header = get_header(SAMReader(IOBuffer(sam_data)));

julia> header.comments[1] |> String
"Æxámplè cøntênt with UTF-8"

julia> header.header.sort_order == SortOrders.unsorted
true
```
"""
struct SAMHeader
    # The entire data of the header, in SAM format
    data::ImmutableMemoryView{UInt8}
    header::Union{Nothing, HeaderHD}
    seqs::Vector{HeaderSQ}
    comments::Vector{ImmutableMemoryView{UInt8}}
    read_groups::Vector{HeaderRG}
    programs::Vector{HeaderPG}
end

Base.propertynames(::SAMHeader) = (:header, :seqs, :comments, :read_groups, :programs)

MemoryViews.MemoryView(x::SAMHeader) = x.data

function Base.get(h::AbstractHeader, tag, default)
    tag = convert(AuxTag, tag)::AuxTag
    for (k, v) in h
        k == tag && return v
    end
    return default
end

function Base.iterate(h::AbstractHeader, state...)
    it = @inline split_each(getfield(h, :fields), UInt8('\t'))
    (f, s) = @something iterate(it, state...) return nothing
    length(f) < 4 && throw(ArgumentError("Malformed AbstractHeader field"))
    tag = AuxTag(@inbounds(f[1]), @inbounds(f[2]))
    val = @inbounds f[4:end]
    return (tag => val, s)
end

Base.length(h::AbstractHeader) = getfield(h, :len)

"""
    try_parse_header(data::ImmutableMemoryView{UInt8})::Union{SAMHeader, SAMHeaderError}

Parse `data` as a SAM header, or return a `SAMHeaderError` with information about how
the parsing failed.
"""
function try_parse_header end

@eval function try_parse_header(data::ImmutableMemoryView{UInt8})
    line_number = 1
    mark = 0
    content_mark = 0
    field = data[1:0]
    field_mem = field
    tags = UInt16[]
    fields = HeaderField[]
    tag = 0x0000
    hd = nothing
    sqs = HeaderSQ[]
    rgs = HeaderRG[]
    pgs = HeaderPG[]
    cos = ImmutableMemoryView{UInt8}[]

    $(Base.remove_linenums!(Automa.generate_init_code(CTX, HEADER_MACHINE)))
    $(Base.remove_linenums!(Automa.generate_exec_code(CTX, HEADER_MACHINE, HEADER_ACTIONS)))

    if cs < 0
        col = p - something(findprev(==(0x0a), data, p), 0)
        return SAMHeaderError(line_number, (col, byte))
    end

    return SAMHeader(data, hd, sqs, cos, rgs, pgs)
end

function has_duplicate_tags(tags::Vector{UInt16})::Bool
    length(tags) < 2 && return false
    if length(tags) < 6
        @inbounds for i in 1:(length(tags) - 1)
            for j in (i + 1):length(tags)
                tags[i] == tags[j] && return true
            end
        end
        return false
    end
    sort!(tags; alg = Base.Sort.QuickSortAlg())
    @inbounds for i in 1:(length(tags) - 1)
        tags[i] == tags[i + 1] && return true
    end
    return false
end

function try_parse_hd(
        fields::Vector{HeaderField},
        fieldmem::ImmutableMemoryView{UInt8}
    )::Union{SAMHeaderError, HeaderHD}
    major = nothing
    minor = UInt32(0)
    sort_order = SortOrders.unknown

    for field in fields
        if field.tag == str_as_u16("VN")
            major_minor = parse_vn_field(field.data)
            if major_minor === nothing
                return SAMHeaderError(1, SAMHeaderErrors.invalid_vn_field)
            end
            (major, minor) = major_minor
        elseif field.tag == str_as_u16("SO")
            sort_order = parse_mem(SortOrder, field.data)
            if sort_order isa SAMHeaderErrorKind
                return SAMHeaderError(1, sort_order)
            end
        end
    end
    major === nothing && return SAMHeaderError(1, SAMHeaderErrors.no_vn_field)
    return HeaderHD(fieldmem, length(fields), major, minor, sort_order)
end

function parse_vn_field(mem::ImmutableMemoryView{UInt8})::Union{Nothing, Tuple{UInt32, UInt32}}
    p = @something findfirst(==(UInt8('.')), mem) return nothing
    major = @something tryparse_base10(UInt32, @inbounds mem[1:(p - 1)]) return nothing
    minor = @something tryparse_base10(UInt32, @inbounds mem[(p + 1):end])  return nothing
    return (major, minor)
end

function try_parse_sq(
        line_number::Int,
        fields::Vector{HeaderField},
        field_mem::ImmutableMemoryView{UInt8},
    )::Union{HeaderSQ, SAMHeaderError}
    sn = nothing
    ln = nothing
    topology = Topologies.linear

    for field_index in eachindex(fields)
        field = @inbounds fields[field_index]
        if field.tag == str_as_u16("SN")
            sn = field.data
        elseif field.tag == str_as_u16("LN")
            ln = tryparse_base10(UInt32, field.data)
            if ln === nothing
                return SAMHeaderError(line_number, SAMHeaderErrors.invalid_ln_field)
            end
        elseif field.tag == str_as_u16("TP")
            topology = if field.data == b"linear"
                Topologies.linear
            elseif field.data == b"circular"
                Topologies.circular
            else
                return SAMHeaderError(line_number, SAMHeaderErrors.invalid_topology)
            end
        end
    end
    ln === nothing && return SAMHeaderError(line_number, SAMHeaderErrors.no_ln_field)
    sn === nothing && return SAMHeaderError(line_number, SAMHeaderErrors.no_sn_field)
    return HeaderSQ(field_mem, sn, length(fields), ln, topology)
end

function try_parse_rg_pg(
        line_number::Int,
        fields::Vector{HeaderField},
    )::Union{SAMHeaderError, ImmutableMemoryView{UInt8}}
    id = nothing
    for field in fields
        if field.tag == str_as_u16("ID")
            id = field.data
        end
    end

    id === nothing && return SAMHeaderError(line_number, SAMHeaderErrors.no_id_field)
    return id
end
