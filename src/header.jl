# Header validation policy: Strictly validate what XAM2 interprets, and store everything
# else raw. Validated:
# * Line structure: `@` and a known record type (HD, SQ, RG, PG, CO). Except for CO,
#   followed by tab-separated `TAG:VALUE` fields with unique tags and non-empty values.
# * HD: Must be the first line. VN present and `major.minor`. SO, if present, is known.
# * SQ: SN present, a valid reference name, and unique across lines. LN present and in
#   1:typemax(Int32). TP, if present, is known.
# * RG and PG: ID present.
# Not validated: Other tags and their values, uniqueness of RG and PG IDs, PG chains (PP).
# The BAM reader additionally checks that the @SQ lines, if any, agree with the binary
# reference list.

@noinline function throw_property_error(x, sym)
    s = "Type $(typeof(x)) has no such property: $(sym)"
    error(s)
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
    `unknown`, `unsorted`, `queryname`, `coordinate`
    """
    @enum SortOrder::UInt8 unknown unsorted queryname coordinate
    export SortOrder, unknown, unsorted, queryname, coordinate
end

using .SortOrders: SortOrders, SortOrder

function parse_sort_order(mem::ImmutableMemoryView{UInt8})::Union{SortOrder, Nothing}
    return if mem == b"unknown"
        SortOrders.unknown
    elseif mem == b"unsorted"
        SortOrders.unsorted
    elseif mem == b"queryname"
        SortOrders.queryname
    elseif mem == b"coordinate"
        SortOrders.coordinate
    else
        nothing
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

const HeaderString = StringView{ImmutableMemoryView{UInt8}}

"""
    AbstractHeader <: AbstractDict{AuxTag, StringView{ImmutableMemoryView{UInt8}}}

A SAM header line, e.g. `@SQ\\tSN:sequence1\\tLN:2015`.
This is an immutable `AbstractDict` where the keys are `XAMAuxData.AuxTag`, and
the values are the corresponding values as `AbstractString`s, which alias the header.

The content of the line, not including the first `@`-prefixed tag and the first
tab, can be obtained with `MemoryView(::AbstractHeader)`.

See also: [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> sam_header = get_header(SAMReader(IOBuffer(sam_data)));

julia> sq_header = sam_header.seqs[1]
XAM2.HeaderSQ with 2 entries:
  AuxTag("SN") => "CHROMOSOME_I"
  AuxTag("LN") => "1009800"

julia> sq_header["LN"]
"1009800"

julia> MemoryView(sq_header) |> String
"SN:CHROMOSOME_I\\tLN:1009800"
```
"""
abstract type AbstractHeader <: AbstractDict{AuxTag, HeaderString} end

function Base.getproperty(x::AbstractHeader, symbol::Symbol)
    if symbol ∈ propertynames(x)
        return getfield(x, symbol)
    end
    throw_property_error(x, symbol)
end

MemoryViews.MemoryView(x::AbstractHeader) = getfield(x, :fields)

# Invariant of every `AbstractHeader`: `fields` is one or more tab-separated fields of
# the form `TAG:VALUE`, where `TAG` is two bytes and `VALUE` is non-empty, and `len` is
# the number of fields. This is checked by `try_parse_header`.

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

julia> hd_header["VN"]
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
* `sn::AbstractString`: Name of sequence
* `ln::Int`: Length of the sequence, in `1:typemax(Int32)`
* `topology::Topology`. Defaults to `Topologies.linear`

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> sq_header = get_header(SAMReader(IOBuffer(sam_data))).seqs[1];

julia> sq_header.topology == Topologies.linear
true

julia> sq_header.ln
1009800

julia> sq_header.sn
"CHROMOSOME_I"
```
"""
struct HeaderSQ <: AbstractHeader
    # All content of the header are stored in contiguous memory
    fields::ImmutableMemoryView{UInt8}
    sn::HeaderString # a view into `fields`
    len::Int # number of fields
    ln::Int
    topology::Topology # linear by default
end

Base.propertynames(::HeaderSQ) = (:sn, :ln, :topology)

"""
    HeaderRG <: AbstractHeader

A struct representing a `@RG` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `id::AbstractString`: Read group identifier

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)
"""
struct HeaderRG <: AbstractHeader
    # All content of the header are stored in contiguous memory
    fields::ImmutableMemoryView{UInt8}
    id::HeaderString # a view into `fields`
    len::Int # number of fields
end

Base.propertynames(::HeaderRG) = (:id,)

"""
    HeaderPG <: AbstractHeader

A struct representing a `@PG` header.
Beside the `AbstractHeader` interface, this type contains the following public properties.
More may be added in the future:
* `id::AbstractString`: Program record identifier

See also: [`AbstractHeader`](@ref), [`SAMHeader`](@ref)

# Examples
```jldoctest
julia> pg_header = get_header(SAMReader(IOBuffer(sam_data))).programs[1];

julia> pg_header.id
"samtools"

julia> MemoryView(pg_header) |> String
"ID:samtools\\tPN:samtools\\tVN:1.23\\tCL:samtools view -h ex1.bam"
```
"""
struct HeaderPG <: AbstractHeader
    # Same fields as HeaderRG
    fields::ImmutableMemoryView{UInt8}
    id::HeaderString
    len::Int
end

Base.propertynames(::HeaderPG) = (:id,)

"""
    SAMHeader

Object representing a SAM or BAM header, i.e. all the lines in a SAM
file starting with `@`.

Obtain a `SAMHeader` from a [`SAMReader`](@ref) or [`BAMReader`](@ref)
using [`get_header`](@ref), or parse one with [`XAM2.try_parse_header`](@ref).

Obtain the header in SAM format with `MemoryView(::SAMHeader)`.

The parts of the header that XAM2 interprets are validated when parsing,
e.g. the `SN` and `LN` fields of `@SQ` lines, and `SN` must be unique.
Other fields are stored as-is, without validation.

A `SAMHeader` has the public properties listed below. More may be added in
future releases. Within each `Vector` of `AbstractHeader`, the elements
appear in the same order as in the source SAM/BAM file.
* `header::Union{Nothing, HeaderHD}`
* `seqs::Vector{HeaderSQ}`
* `comments::Vector{<:AbstractString}`: This represents all the `@CO`
  lines, without the leading `"@CO\\t"`.
* `read_groups::Vector{HeaderRG}`
* `programs::Vector{HeaderPG}`

See also: [`AbstractHeader`](@ref), [`HeaderHD`](@ref)

# Examples
```jldoctest
julia> header = get_header(SAMReader(IOBuffer(sam_data)));

julia> header.comments[1]
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
    comments::Vector{HeaderString}
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
    # Invariant: Every field is TAG:VALUE, with a two-byte tag
    tag = AuxTag(@inbounds(f[1]), @inbounds(f[2]))
    val = StringView(@inbounds f[4:end])
    return (tag => val, s)
end

Base.length(h::AbstractHeader) = getfield(h, :len)

"""
    try_parse_header(data::ImmutableMemoryView{UInt8})::Union{SAMHeader, XAMError}

Parse `data` as a SAM header, or return an `XAMError` with the kind of error, and
the line and column where it occurred. Lines may end with `\\n` or `\\r\\n`, and the last line need not
end with a newline. The returned header aliases `data`.
"""
function try_parse_header(data::ImmutableMemoryView{UInt8})::Union{SAMHeader, XAMError}
    state = HeaderParseState()
    for (line_number, line) in enumerate(each_header_line(data))
        err = parse_header_line!(state, line, line_number)
        if err !== nothing
            (column, kind) = err
            return XAMError(kind; line = line_number, column)
        end
    end
    return SAMHeader(data, state.hd, state.sqs, state.cos, state.rgs, state.pgs)
end

# Iterator of the lines of SAM header text, without the line terminators,
# which are `\n` or `\r\n`. The last line need not end with a newline.
function each_header_line(data::ImmutableMemoryView{UInt8})
    # A final newline ends the last line, rather than starting an empty line.
    # But if the newline is the only byte, it ends an empty line, which is invalid.
    if length(data) > 1 && @inbounds(data[end]) == UInt8('\n')
        data = @inbounds data[1:(end - 1)]
    end
    return Iterators.map(strip_carriage_return, split_each(data, UInt8('\n')))
end

function strip_carriage_return(line::ImmutableMemoryView{UInt8})::ImmutableMemoryView{UInt8}
    return if !isempty(line) && @inbounds(line[end]) == UInt8('\r')
        @inbounds line[1:(end - 1)]
    else
        line
    end
end

struct HeaderField
    tag::UInt16 # the two-byte tag stored in little-endian format
    column::Int # column of the start of the field in the line
    value::ImmutableMemoryView{UInt8}
end

mutable struct HeaderParseState
    hd::Union{Nothing, HeaderHD}
    const sqs::Vector{HeaderSQ}
    const rgs::Vector{HeaderRG}
    const pgs::Vector{HeaderPG}
    const cos::Vector{HeaderString}
    # SN fields seen so far, to check they are unique
    const sns::Set{HeaderString}
    # Fields of the current line
    const fields::Vector{HeaderField}
    # Bit `i` is set if two-byte tag `i - 1` was seen in the current line. Cleared
    # after each line, by unsetting the tags in `fields`.
    const seen_tags::BitVector
end

function HeaderParseState()
    return HeaderParseState(
        nothing, HeaderSQ[], HeaderRG[], HeaderPG[], HeaderString[],
        Set{HeaderString}(), HeaderField[], falses(2^16),
    )
end

str_as_u16(s::String) = (codeunit(s, 1) % UInt16) | ((codeunit(s, 2) % UInt16) << 8)

# Errors are returned as (column, kind)
const HeaderLineError = Tuple{Int, XAMErrorKind}

function parse_header_line!(
        state::HeaderParseState,
        line::ImmutableMemoryView{UInt8},
        line_number::Int,
    )::Union{Nothing, HeaderLineError}
    # The record type, and if there is content, a tab after it
    if length(line) < 3 || @inbounds(line[1]) != UInt8('@') ||
            (length(line) > 3 && @inbounds(line[4]) != UInt8('\t'))
        return (1, XAMErrors.header_invalid_line)
    end
    record_type = (@inbounds(line[2]) % UInt16) | ((@inbounds(line[3]) % UInt16) << 8)
    # Content after the record type and the tab
    content = @inbounds line[min(5, length(line) + 1):end]

    if record_type == str_as_u16("CO")
        # Comments are arbitrary text
        push!(state.cos, StringView(content))
        return nothing
    elseif record_type ∉ (str_as_u16("HD"), str_as_u16("SQ"), str_as_u16("RG"), str_as_u16("PG"))
        return (2, XAMErrors.header_unknown_record_type)
    end

    err = split_header_fields!(state, line)
    err === nothing || return err
    fields = state.fields

    return if record_type == str_as_u16("HD")
        line_number == 1 || return (1, XAMErrors.header_hd_not_first)
        hd = parse_hd(fields, content)
        hd isa HeaderLineError && return hd
        state.hd = hd
        nothing
    elseif record_type == str_as_u16("SQ")
        sq = parse_sq!(state.sns, fields, content)
        sq isa HeaderLineError && return sq
        push!(state.sqs, sq)
        nothing
    else
        id = parse_id(fields)
        id isa HeaderLineError && return id
        if record_type == str_as_u16("RG")
            push!(state.rgs, HeaderRG(content, id, length(fields)))
        else
            push!(state.pgs, HeaderPG(content, id, length(fields)))
        end
        nothing
    end
end

# Split the TAG:VALUE fields of a non-comment line into `state.fields`,
# checking their syntax, and that tags are unique.
function split_header_fields!(
        state::HeaderParseState,
        line::ImmutableMemoryView{UInt8},
    )::Union{Nothing, HeaderLineError}
    fields = state.fields
    seen = state.seen_tags
    # Clear the tags of the previous line
    for field in fields
        @inbounds seen[(field.tag % Int) + 1] = false
    end
    empty!(fields)
    # Index of the tab preceding the current field
    tab = 4
    while tab ≤ length(line)
        column = tab + 1
        next_tab = something(findnext(==(UInt8('\t')), line, column), length(line) + 1)
        field = @inbounds line[column:(next_tab - 1)]
        is_valid_field(field) || return (column, XAMErrors.header_invalid_field)
        tag = (@inbounds(field[1]) % UInt16) | ((@inbounds(field[2]) % UInt16) << 8)
        @inbounds(seen[(tag % Int) + 1]) && return (column, XAMErrors.header_duplicate_tag)
        @inbounds seen[(tag % Int) + 1] = true
        push!(fields, HeaderField(tag, column, @inbounds field[4:end]))
        tab = next_tab
    end
    return nothing
end

function is_valid_field(field::ImmutableMemoryView{UInt8})::Bool
    length(field) < 4 && return false
    (a, b, c) = @inbounds (field[1], field[2], field[3])
    is_letter(x) = in(x, UInt8('A'):UInt8('Z')) | in(x, UInt8('a'):UInt8('z'))
    return is_letter(a) & (is_letter(b) | in(b, UInt8('0'):UInt8('9'))) & (c == UInt8(':'))
end

function parse_hd(
        fields::Vector{HeaderField},
        content::ImmutableMemoryView{UInt8},
    )::Union{HeaderHD, HeaderLineError}
    version = nothing
    sort_order = SortOrders.unknown
    for field in fields
        if field.tag == str_as_u16("VN")
            version = parse_vn_field(field.value)
            version === nothing && return (field.column, XAMErrors.header_invalid_vn)
        elseif field.tag == str_as_u16("SO")
            sort_order = @something parse_sort_order(field.value) begin
                return (field.column, XAMErrors.header_invalid_sort_order)
            end
        end
    end
    version === nothing && return (1, XAMErrors.header_missing_vn)
    (major, minor) = version
    return HeaderHD(content, length(fields), major, minor, sort_order)
end

function parse_vn_field(mem::ImmutableMemoryView{UInt8})::Union{Nothing, Tuple{UInt32, UInt32}}
    p = @something findfirst(==(UInt8('.')), mem) return nothing
    major = @something tryparse_base10(UInt32, @inbounds mem[1:(p - 1)]) return nothing
    minor = @something tryparse_base10(UInt32, @inbounds mem[(p + 1):end]) return nothing
    return (major, minor)
end

# Parse an @SQ line, and add its name to `sns`, the names seen in previous lines.
function parse_sq!(
        sns::Set{HeaderString},
        fields::Vector{HeaderField},
        content::ImmutableMemoryView{UInt8},
    )::Union{HeaderSQ, HeaderLineError}
    sn = nothing
    sn_column = 0
    ln = nothing
    topology = Topologies.linear

    for field in fields
        if field.tag == str_as_u16("SN")
            is_valid_reference_name(field.value) ||
                return (field.column, XAMErrors.header_invalid_sn)
            sn = StringView(field.value)
            sn_column = field.column
        elseif field.tag == str_as_u16("LN")
            ln = tryparse_base10(Int32, field.value)
            # Note: This also rejects a leading '-', as that would give a negative number
            (ln === nothing || ln < 1) && return (field.column, XAMErrors.header_invalid_ln)
        elseif field.tag == str_as_u16("TP")
            topology = if field.value == b"linear"
                Topologies.linear
            elseif field.value == b"circular"
                Topologies.circular
            else
                return (field.column, XAMErrors.header_invalid_topology)
            end
        end
    end
    sn === nothing && return (1, XAMErrors.header_missing_sn)
    ln === nothing && return (1, XAMErrors.header_missing_ln)
    n_before = length(sns)
    push!(sns, sn)
    length(sns) == n_before && return (sn_column, XAMErrors.header_duplicate_sn)
    return HeaderSQ(content, sn, length(fields), ln % Int, topology)
end

function parse_id(fields::Vector{HeaderField})::Union{HeaderString, HeaderLineError}
    for field in fields
        field.tag == str_as_u16("ID") && return StringView(field.value)
    end
    return (1, XAMErrors.header_missing_id)
end

# Bytes allowed in reference names, see the SAM specification, section 1.2.1:
# [0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*
const REFERENCE_NAME_BYTES = let
    bits = falses(256)
    for byte in b"!#$%&+./:;?@^_|~-*="
        bits[byte + 1] = true
    end
    for range in Any['A':'Z', 'a':'z', '0':'9'], char in range
        bits[UInt8(char) + 1] = true
    end
    Tuple(bits)
end

function is_valid_reference_name(name::ImmutableMemoryView{UInt8})::Bool
    isempty(name) && return false
    first_byte = @inbounds name[1]
    (first_byte == UInt8('*')) | (first_byte == UInt8('=')) && return false
    valid = true
    for byte in name
        valid &= @inbounds REFERENCE_NAME_BYTES[byte + 1]
    end
    return valid
end
