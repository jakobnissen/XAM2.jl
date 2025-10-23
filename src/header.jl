const HEADER_MACHINE = let
    lf = onenter!(re"\n", :newline)
    newline = opt(re"\r") * lf
    header_field = onexit!(onenter!(re"[ -~]+", :mark), :field)
    tag = onexit!(re"[A-Za-z][A-Za-z0-9]", :tag)
    header_segment = onexit!(tag * re":" * header_field, :segment)
    header_content = header_segment * rep(re"\t" * header_segment)
    onexit!(onenter!(header_content, :mark_content), :header_content)

    header_start = onenter!(re"@", :header)
    # TODO: Should be able to contain \r, but I disallow it here for Automa amgiguitity reasons
    hd = onexit!(header_start * re"HD\t" * header_content, [:check_unique, :hd])
    co = onexit!(header_start * onenter!(re"CO\t[^\n\r]*", :enter_co), :co)
    sq = onexit!(header_start * re"SQ\t" * header_content, [:check_unique, :sq])
    rg = onexit!(header_start * re"RG\t" * header_content, [:check_unique, :rg])
    pg = onexit!(header_start * re"PG\t" * header_content, [:check_unique, :pg])

    header = opt(hd) * rep(newline * (co | sq | rg | pg))
    Automa.compile(header * rep(newline))
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
    :check_unique => :(check_duplicate_tags(line_number, tags)),
    :hd => :(hd = parse_hd(fields, field_mem)),
    :enter_co => :(content_mark = p + 3),
    :co => :(push!(cos, @inbounds(data[content_mark:(p - 1)]))),
    :sq => :(push!(sqs, parse_sq(line_number, fields, field_mem))),
    :rg => quote
        id = parse_rg_pg(line_number, fields)
        push!(rgs, HeaderRG(field_mem, id, length(fields)))
    end,
    :pg => quote
        id = parse_rg_pg(line_number, fields)
        push!(pgs, HeaderPG(field_mem, id, length(fields)))
    end,
)

struct HeaderError <: Exception
    line::Int
    str::String
end

struct HeaderField
    tag::UInt16
    data::ImmutableMemoryView{UInt8}
end

function str_as_u16(s)
    ncodeunits(s) == 2 || error()
    return (@inbounds(codeunit(s, 1)) % UInt16) | ((@inbounds(codeunit(s, 2)) % UInt16) << 8)
end

module SortOrders
    @enum SortOrder::UInt8 unknown unsorted queryname coordinate
    export SortOrder, unknown, unsorted, queryname, coordinate
end

using .SortOrders: SortOrders, SortOrder

function parse_mem(::Type{SortOrder}, mem::ImmutableMemoryView{UInt8})
    return if mem == b"unknown"
        SortOrders.unknown
    elseif mem == b"unsorted"
        SortOrders.unsorted
    elseif mem == b"queryname"
        SortOrders.queryname
    elseif mem == b"coordinate"
        SortOrders.coordinate
    else
        # HD line must be line 1
        throw(HeaderError(1, "In HD line, SO field is not \"unknown\", \"unsorted\", \"queryname\" or \"coordinate\""))
    end
end

module Topologies
    @enum Topology::UInt8 linear circular
    export Topology, linear, circular
end

using .Topologies: Topologies, Topology

abstract type AbstractHeader <: AbstractDict{AuxTag, ImmutableMemoryView{UInt8}} end

struct HeaderHD <: AbstractHeader
    fields::ImmutableMemoryView{UInt8}
    len::Int
    major::UInt32
    minor::UInt32
    sort_order::SortOrder
end

struct HeaderSQ <: AbstractHeader
    fields::ImmutableMemoryView{UInt8}
    sn::ImmutableMemoryView{UInt8}
    len::Int
    ln::UInt32
    topology::Topology # linear by default
end

struct HeaderRG <: AbstractHeader
    fields::ImmutableMemoryView{UInt8}
    id::ImmutableMemoryView{UInt8}
    len::Int
end

struct HeaderPG <: AbstractHeader
    fields::ImmutableMemoryView{UInt8}
    id::ImmutableMemoryView{UInt8}
    len::Int
end

# TODO: Store the overall view somewhere in this struct, then impl MemoryView(::SAMHeader) as immut view
struct SAMHeader
    header::Union{Nothing, HeaderHD}
    seqs::Vector{HeaderSQ}
    comments::Vector{ImmutableMemoryView{UInt8}}
    read_groups::Vector{HeaderRG}
    programs::Vector{HeaderPG}
end

function Base.get(h::AbstractHeader, tag::AuxTag, default)
    for (k, v) in h
        k == tag && return v
    end
    return default
end

function Base.iterate(h::AbstractHeader, state...)
    it = @inline split_each(h.fields, UInt8('\t'))
    (f, s) = @something iterate(it, state...) return nothing
    length(f) < 4 && throw(ArgumentError("Malformed AbstractHeader field"))
    tag = AuxTag(@inbounds(f[1]), @inbounds(f[2]))
    val = @inbounds f[4:end]
    return (tag => val, s)
end

Base.length(h::AbstractHeader) = h.len

@eval function parse_header(data::ImmutableMemoryView{UInt8})
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
        throw_automa_error(byte, line_number, col)
    end

    return SAMHeader(hd, sqs, cos, rgs, pgs)
end

@noinline function throw_automa_error(byte::UInt8, line_number::Int, col::Int)
    msg = string(
        "Cannot parse SAM/BAM header, got unexpected byte ",
        repr(byte),
        " at line ",
        line_number,
        " column ",
        col,
        "."
    )
    throw(HeaderError(line_number, msg))
end

function check_duplicate_tags(line_number::Int, tags::Vector{UInt16})
    length(tags) < 2 && return nothing
    if length(tags) < 6
        @inbounds for i in 1:(length(tags) - 1)
            for j in (i + 1):length(tags)
                tags[i] == tags[j] && @goto throw_duplicate
            end
        end
        return
    end
    sort!(tags; alg = Base.Sort.QuickSortAlg())
    @inbounds for i in 1:(length(tags) - 1)
        tags[i] == tags[i + 1] && @goto throw_duplicate
    end
    return nothing
    @label throw_duplicate
    throw(HeaderError(line_number, "Duplicate tag found"))
end

function parse_hd(fields::Vector{HeaderField}, fieldmem::ImmutableMemoryView{UInt8})::HeaderHD
    major = nothing
    minor = UInt32(0)
    sort_order = SortOrders.unknown

    for field in fields
        if field.tag == str_as_u16("VN")
            major_minor = parse_vn_field(field.data)
            major_minor === nothing && throw(
                HeaderError(
                    "In HD line, VN field is not formatted [0-9]+\\.[0-9]+, or either version number overflowed."
                )
            )
            (major, minor) = major_minor
        elseif field.tag == str_as_u16("SO")
            sort_order = parse_mem(SortOrder, field.data)
        end
    end
    major === nothing && throw(HeaderError("In HD line, did not find mandatory VN field"))
    return HeaderHD(fieldmem, length(fields), major, minor, sort_order)
end

function parse_vn_field(mem::ImmutableMemoryView{UInt8})::Union{Nothing, Tuple{UInt32, UInt32}}
    p = @something findfirst(==(UInt8('.')), mem) return nothing
    major = @something tryparse_base10(UInt32, @inbounds mem[1:(p - 1)]) return nothing
    minor = @something tryparse_base10(UInt32, @inbounds mem[(p + 1):end])  return nothing
    return (major, minor)
end

function parse_sq(
        line_number::Int,
        fields::Vector{HeaderField},
        field_mem::ImmutableMemoryView{UInt8},
    )::HeaderSQ
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
                throw(HeaderError(line_number, "Could not parse LN field in SQ header as 32-bit integer"))
            end
        elseif field.tag == str_as_u16("TP")
            topology = if field.data == b"linear"
                Topologies.linear
            elseif field.data == b"circular"
                Topologies.circular
            else
                throw(HeaderError(line_number, "TP field in SQ line must be \"linear\" or \"circular\""))
            end
        end
    end
    ln === nothing && throw(HeaderError(line_number, "SQ line does not have mandatory LN field"))
    sn === nothing && throw(HeaderError(line_number, "SQ line does not have mandatory SN field"))
    return HeaderSQ(field_mem, sn, length(fields), ln, topology)
end

function parse_rg_pg(
        line_number::Int,
        fields::Vector{HeaderField},
    )::ImmutableMemoryView{UInt8}
    id = nothing
    for field in fields
        if field.tag == str_as_u16("ID")
            id = field.data
        end
    end
    id === nothing && throw(HeaderError(line_number, "RG or PG header line must have ID field"))
    return id
end
