const TAG_HD = AuxTag("HD")
const TAG_SQ = AuxTag("SQ")
const TAG_RG = AuxTag("RG")
const TAG_PG = AuxTag("PG")
const TAG_CO = AuxTag("CO")

struct HeaderError <: Exception
    line::Int
    str::String
end

struct HeaderField
    tag::AuxTag
    data::ImmutableMemoryView{UInt8}
end

module SortOrders
    @enum SortOrder::UInt16 unknown unsorted queryname coordinate
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

# TODO: Do not store a vector of HeaderFIeld, jsut store the raw
# data
struct HeaderHD <: AbstractHeader
    fields::Vector{HeaderField}
    major::UInt32
    minor::UInt32
    sort_order::SortOrder
end

struct HeaderSQ <: AbstractHeader
    fields::Vector{HeaderField}
    sn::ImmutableMemoryView{UInt8}
    ln::UInt32
    topology::Topology # linear by default
end

struct HeaderRG <: AbstractHeader
    fields::Vector{HeaderField}
    id::ImmutableMemoryView{UInt8}
end

struct HeaderPG <: AbstractHeader
    fields::Vector{HeaderField}
    id::ImmutableMemoryView{UInt8}
end

struct HeaderOther <: AbstractHeader
    tag::AuxTag
    fields::Vector{HeaderField}
end

struct SAMHeader
    header::Union{Nothing, HeaderHD}
    seqs::Vector{HeaderSQ}
    comments::Vector{ImmutableMemoryView{UInt8}}
    read_groups::Vector{HeaderRG}
    programs::Vector{HeaderPG}
    other::Vector{HeaderOther}
end

function Base.get(h::AbstractHeader, tag::AuxTag, default)
    for i in h.fields
        i.tag == tag && return i.data
    end
    return default
end

function Base.iterate(h::AbstractHeader, state...)
    (f, s) = @something iterate(h.fields, state...) return nothing
    return (f.tag => f.data, s)
end

Base.length(h::AbstractHeader) = length(h.fields)

#=
version(::SAMHeader)::Union{Nothing, Tuple{UInt32, UInt32}}
sequences(::SAMHeader)::Vector{HeaderSQ}

Specific headers are AbstractDict, and also have some
public fields.
Fields are also part of the dict (for forward compatibility)






=#

function read_header_if_present!(lines::BufferIO.LineViewIterator{CursorReader})::Union{HeaderHD, Nothing}
    it = iterate(lines)
    isnothing(it) && return nothing
    (memory, state) = it
    length(memory) < 3 && return nothing
    @inbounds(memory[1:3]) == b"@HD" || return nothing

    # Minimum header: "@HD\tVN:1.0" (10 bytes)
    length(memory) < 10 && @goto throw_header
    memory[4] == UInt8('\t') || @goto throw_header

    fields = HeaderField[]
    tags = AuxTag[]
    rest = memory[5:end]
    parse_tags!(1, tags, fields, rest)
    header = parse_header_hd(fields)
    iterate(lines, state) # advance to next line

    # Now, check the version is parseable
    # we currently have version 1.6
    header.major != 1 && throw(HeaderError(1, "XAM.jl only support SAM version 1, got different version in HD line"))
    return header

    @label throw_header
    throw(HeaderError(1, "Malformed HD header line"))
end

function parse_sam_header(data::ImmutableMemoryView{UInt8})
    seqs = HeaderSQ[]
    read_groups = HeaderRG[]
    comments = ImmutableMemoryView{UInt8}[]
    programs = HeaderPG[]
    other = HeaderOther[]
    tags = AuxTag[]
    lines = line_views(CursorReader(data))
    header = read_header_if_present!(lines)

    for (line_index, line) in enumerate(lines)
        length(line) < 3 && throw(HeaderError(line_index, "Malformed header line, less than 3 bytes in length"))
        @inbounds line[1] == UInt8('@') || throw(HeaderError(line_index, "Malformed header line; does not start with @"))

        tag = try_auxtag(@inbounds(line[2]), @inbounds(line[3]))
        tag === nothing && throw(HeaderError(line_index, "Could not parse header tag, must be [A-Za-z][A-Za-z0-9]"))

        if tag == TAG_CO
            rest = @inbounds line[3:end]
            push!(comments, rest)
            continue
        end

        rest = if length(line) == 3
            @inbounds line[4:end]
        else
            if @inbounds(line[4]) != UInt8('\t')
                throw(HeaderError(line_index, "Header line's 4th byte must be a tab"))
            end
            @inbounds line[5:end]
        end

        fields = HeaderField[]
        parse_tags!(line_index, tags, fields, rest)

        if tag == TAG_HD
            throw(HeaderError(line_index, "Found @HD at another line than the first line"))
        elseif tag == TAG_SQ
            seq = parse_header_sq(line_index, fields)
            push!(seqs, seq)
        elseif tag == TAG_RG
            id = parse_header_rg_pg(line_index, fields)
            push!(read_groups, HeaderRG(fields, id))
        elseif tag == TAG_PG
            id = parse_header_rg_pg(line_index, fields)
            push!(programs, HeaderPG(fields, id))
        else
            push!(other, HeaderOther(tag, fields))
        end
    end

    if !isempty(read_groups) || !isempty(programs)
        ids = Set{ImmutableMemoryView{UInt8}}()
        for i in read_groups
            in!(i.id, ids) && throw(HeaderError("ID fields in RG lines are not unique"))
        end
        empty!(ids)
        for i in programs
            in!(i.id, ids) && throw(HeaderError("ID fields in PG lines are not unique"))
        end
    end

    return SAMHeader(header, seqs, comments, read_groups, programs, other)
end

function parse_tags!(
        line_number::Int,
        tags::Vector{AuxTag},
        v::Vector{HeaderField},
        mem::ImmutableMemoryView{UInt8}
    )
    isempty(mem) && return nothing
    empty!(tags)
    for field in split_each(mem, UInt8('\t'))
        length(field) < 3 && @goto throw_bad_tag
        tag = try_auxtag(@inbounds(field[1]), @inbounds(field[2]))
        tag === nothing && @goto throw_bad_tag
        push!(tags, tag)
        @inbounds field[3] == UInt8(':') || @goto throw_bad_tag
        field_mem = @inbounds split_at(field, 4)[2]
        push!(v, HeaderField(tag, field_mem))
    end
    length(tags) < 2 && return nothing
    check_duplicate_tags(line_number, tags)
    return nothing
    @label throw_bad_tag
    throw(HeaderError(line_number, "Field in header line too short, must start with \"[A-Za-z][A-Za-z0-9]:\""))
end

function check_duplicate_tags(line_number::Int, tags::Vector{AuxTag})
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

function parse_header_hd(fields::Vector{HeaderField})::HeaderHD
    major = nothing
    minor = UInt32(0)
    sort_order = SortOrders.unknown

    for field_index in eachindex(fields)
        field = @inbounds fields[field_index]
        if field.tag == try_auxtag("VN")
            major_minor = parse_vn_field(field.data)
            major_minor === nothing && throw(
                HeaderError(
                    "In HD line, VN field is not formatted [0-9]+\\.[0-9]+, or either version number overflowed."
                )
            )
            (major, minor) = major_minor
        elseif field.tag == try_auxtag("SO")
            sort_order = parse_mem(SortOrder, field.data)
        end
    end
    major === nothing && throw(HeaderError("In HD line, did not find mandatory VN field"))
    return HeaderHD(fields, major, minor, sort_order)
end

function parse_vn_field(mem::ImmutableMemoryView{UInt8})::Union{Nothing, Tuple{UInt32, UInt32}}
    p = @something findfirst(==(UInt8('.')), mem) return nothing
    major = @something tryparse_base10(UInt32, @inbounds mem[1:(p - 1)]) return nothing
    minor = @something tryparse_base10(UInt32, @inbounds mem[(p + 1):end])  return nothing
    return (major, minor)
end

function parse_header_sq(line_number::Int, fields::Vector{HeaderField})::Union{Nothing, HeaderSQ}
    sn = nothing
    ln = nothing
    topology = Topologies.linear

    for field_index in eachindex(fields)
        field = @inbounds fields[field_index]
        if field.tag == try_auxtag("SN")
            sn = field.data
        elseif field.tag == try_auxtag("LN")
            ln = tryparse_base10(UInt32, field.data)
            if ln === nothing
                throw(HeaderError(line_number, "Could not parse LN field in SQ header as 32-bit integer"))
            end
        elseif field.tag == try_auxtag("TP")
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
    return HeaderSQ(fields, sn, ln, topology)
end

function parse_header_rg_pg(
        line_number::Int,
        fields::Vector{HeaderField}
    )::ImmutableMemoryView{UInt8}
    id = nothing
    for field in fields
        if field.tag == AuxTag("ID")
            id = field.data
        end
    end
    id === nothing && throw(HeaderError(line_number, "RG or PG header line must have ID field"))
    return id
end
