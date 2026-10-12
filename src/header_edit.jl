# Functions that return an edited copy of a `SAMHeader`.
# Each function edits the lines of the header, then writes and re-parses the whole header,
# so it takes time proportional to the size of the header. This is fine for the few edits
# a typical tool makes, but building a large header with repeated edits is quadratic.

const RECORD_TYPE_ORDER = map(str_as_u16, ("HD", "SQ", "RG", "PG", "CO"))

# Record type of a line of a valid header, encoded like `str_as_u16`
function record_type(line::ImmutableMemoryView{UInt8})::UInt16
    # Invariant: Lines of a valid header start with `@` and a two-letter record type
    return (@inbounds(line[2]) % UInt16) | ((@inbounds(line[3]) % UInt16) << 8)
end

# Position of a record type in the conventional order of header lines: HD, SQ, RG, PG, CO
record_type_rank(type::UInt16)::Int = something(findfirst(==(type), RECORD_TYPE_ORDER))

header_lines(header::SAMHeader) = collect(ImmutableMemoryView{UInt8}, each_header_line(MemoryView(header)))

# Insert `new_lines` of record type `type` after the last line of the same record type,
# or of a record type that conventionally precedes it.
function insert_lines!(
        lines::Vector{ImmutableMemoryView{UInt8}},
        type::UInt16,
        new_lines::Vector{ImmutableMemoryView{UInt8}},
    )::Vector{ImmutableMemoryView{UInt8}}
    rank = record_type_rank(type)
    i = something(findlast(line -> record_type_rank(record_type(line)) ≤ rank, lines), 0)
    splice!(lines, (i + 1):i, new_lines)
    return lines
end

# Join the lines into a new header
function build_header(lines::Vector{ImmutableMemoryView{UInt8}})::SAMHeader
    data = Memory{UInt8}(undef, sum(line -> length(line) + 1, lines; init = 0))
    mem = MemoryView(data)
    i = 0
    for line in lines
        copyto!(@inbounds(mem[(i + 1):(i + length(line))]), line)
        i += length(line) + 1
        @inbounds mem[i] = UInt8('\n')
    end
    result = try_parse_header(ImmutableMemoryView(data))
    # The edited lines are validated by the callers, so this is a bug in XAM2
    result isa XAMError && error("Edited SAM header is invalid: ", sprint(showerror, result))
    return result
end

function check_field_value(tag::String, value::AbstractString)::Nothing
    isempty(value) && throw(ArgumentError(lazy"Value of $tag field is empty"))
    if any(in(('\t', '\n', '\r', '\0')), value)
        throw(
            ArgumentError(
                lazy"Value of $tag field contains a tab, newline, carriage return or null byte"
            )
        )
    end
    return nothing
end

# Print the `tag => value` fields, skipping values that are `nothing`
function print_fields(io::IO, fields)::Nothing
    for (tag, value) in fields
        value === nothing && continue
        check_field_value(tag, value)
        print(io, '\t', tag, ':', value)
    end
    return nothing
end

"""
    add_program(
        header::SAMHeader, id::AbstractString;
        pn=nothing, cl=nothing, pp=<last program>, ds=nothing, vn=nothing,
    )::SAMHeader

Return a copy of `header` with a new `@PG` line, placed after the existing `@PG` lines.

The keyword arguments are the optional fields of the line: `pn` (program name), `cl`
(command line), `pp` (previous program), `ds` (description) and `vn` (program version).
Their values are `AbstractString`s, or `nothing` to omit the field.

If `id` is the ID of an existing `@PG` line, the suffix `.1`, `.2` etc. is appended to
make it unique, as samtools does. Get the ID of the new line with
`last(new_header.programs).id`.

By default, `pp` is the ID of the last program in the chain of `@PG` lines, i.e. the line
whose ID is not the `PP` of any other `@PG` line. If there are no `@PG` lines, or several
chains, e.g. from merging files, `pp` defaults to `nothing`.

Throw an `ArgumentError` if `pp` is not the ID of a `@PG` line in `header`, or if `id` or
a field value is empty, or contains a tab, newline, carriage return or null byte.

This takes time proportional to the size of `header`. To build a header with many lines,
it is more efficient to write its SAM text, and parse it with [`XAM2.try_parse_header`](@ref).

See also: [`XAM2.add_read_group`](@ref)

# Examples
```jldoctest
julia> header = XAM2.try_parse_header(ImmutableMemoryView(
           "@HD\\tVN:1.6\\n@PG\\tID:tool\\n@CO\\tA comment\\n"
       ));

julia> header = XAM2.add_program(header, "tool"; pn = "tool", vn = "1.2");

julia> print(String(MemoryView(header)))
@HD	VN:1.6
@PG	ID:tool
@PG	ID:tool.1	PN:tool	PP:tool	VN:1.2
@CO	A comment
```
"""
function add_program(
        header::SAMHeader,
        id::AbstractString;
        pn::Union{Nothing, AbstractString} = nothing,
        cl::Union{Nothing, AbstractString} = nothing,
        pp::Union{Nothing, AbstractString} = last_program_id(header),
        ds::Union{Nothing, AbstractString} = nothing,
        vn::Union{Nothing, AbstractString} = nothing,
    )::SAMHeader
    check_field_value("ID", id)
    ids = Set{HeaderString}(pg.id for pg in header.programs)
    if pp !== nothing && pp ∉ ids
        throw(ArgumentError(lazy"PP \"$pp\" is not the ID of a @PG line in the header"))
    end
    io = IOBuffer()
    print(io, "@PG\tID:", unique_program_id(ids, id))
    print_fields(io, ("PN" => pn, "CL" => cl, "PP" => pp, "DS" => ds, "VN" => vn))
    lines = header_lines(header)
    insert_lines!(lines, str_as_u16("PG"), [ImmutableMemoryView(take!(io))])
    return build_header(lines)
end

# The ID of the @PG line that is not the PP of any other @PG line, or `nothing`
# if there are zero or several such IDs.
function last_program_id(header::SAMHeader)::Union{Nothing, HeaderString}
    previous = Set{HeaderString}()
    for pg in header.programs
        pp = get(pg, "PP", nothing)
        pp === nothing || push!(previous, pp)
    end
    result = nothing
    for pg in header.programs
        pg.id ∈ previous && continue
        if result === nothing
            result = pg.id
        elseif result != pg.id
            return nothing
        end
    end
    return result
end

function unique_program_id(ids::Set{HeaderString}, id::AbstractString)::String
    id ∈ ids || return String(id)
    n = 1
    while string(id, '.', n) ∈ ids
        n += 1
    end
    return string(id, '.', n)
end

"""
    add_read_group(
        header::SAMHeader, id::AbstractString;
        bc=nothing, cn=nothing, ds=nothing, dt=nothing, fo=nothing, ks=nothing, lb=nothing,
        pg=nothing, pi=nothing, pl=nothing, pm=nothing, pu=nothing, sm=nothing,
    )::SAMHeader

Return a copy of `header` with a new `@RG` line with the ID `id`, placed after the existing
`@RG` lines.

The keyword arguments are the optional fields of the line, e.g. `sm` (sample) and `lb`
(library). See the SAM specification for their meaning. Their values are `AbstractString`s,
or `nothing` to omit the field. Their values are not validated, beyond the checks below.

Throw an `ArgumentError` if `header` already has a `@RG` line with the ID `id`, or if
`id` or a field value is empty, or contains a tab, newline, carriage return or null byte.

This takes time proportional to the size of `header`. To build a header with many lines,
it is more efficient to write its SAM text, and parse it with [`XAM2.try_parse_header`](@ref).

See also: [`XAM2.add_program`](@ref)

# Examples
```jldoctest
julia> header = XAM2.try_parse_header(ImmutableMemoryView("@SQ\\tSN:chr1\\tLN:100\\n"));

julia> header = XAM2.add_read_group(header, "rg1"; sm = "sample1", pl = "ILLUMINA");

julia> only(header.read_groups)["SM"]
"sample1"
```
"""
function add_read_group(
        header::SAMHeader,
        id::AbstractString;
        bc::Union{Nothing, AbstractString} = nothing,
        cn::Union{Nothing, AbstractString} = nothing,
        ds::Union{Nothing, AbstractString} = nothing,
        dt::Union{Nothing, AbstractString} = nothing,
        fo::Union{Nothing, AbstractString} = nothing,
        ks::Union{Nothing, AbstractString} = nothing,
        lb::Union{Nothing, AbstractString} = nothing,
        pg::Union{Nothing, AbstractString} = nothing,
        pi::Union{Nothing, AbstractString} = nothing,
        pl::Union{Nothing, AbstractString} = nothing,
        pm::Union{Nothing, AbstractString} = nothing,
        pu::Union{Nothing, AbstractString} = nothing,
        sm::Union{Nothing, AbstractString} = nothing,
    )::SAMHeader
    check_field_value("ID", id)
    # Records refer to read groups by ID, so a duplicate can't be renamed like @PG IDs
    if any(rg -> rg.id == id, header.read_groups)
        throw(ArgumentError(lazy"Header already has a @RG line with ID \"$id\""))
    end
    io = IOBuffer()
    print(io, "@RG\tID:", id)
    fields = (
        "BC" => bc, "CN" => cn, "DS" => ds, "DT" => dt, "FO" => fo, "KS" => ks, "LB" => lb,
        "PG" => pg, "PI" => pi, "PL" => pl, "PM" => pm, "PU" => pu, "SM" => sm,
    )
    print_fields(io, fields)
    lines = header_lines(header)
    insert_lines!(lines, str_as_u16("RG"), [ImmutableMemoryView(take!(io))])
    return build_header(lines)
end

function sort_order_name(x::SortOrder)::String
    return if x == SortOrders.unknown
        "unknown"
    elseif x == SortOrders.unsorted
        "unsorted"
    elseif x == SortOrders.queryname
        "queryname"
    else
        @assert x == SortOrders.coordinate
        "coordinate"
    end
end

"""
    set_sort_order(header::SAMHeader, sort_order::SortOrder)::SAMHeader

Return a copy of `header` where the `SO` field of the `@HD` line is `sort_order`.
The `SS` (sub-sorting) and `GO` (grouping) fields are removed, since they describe the
previous order of the records.
If `header` has no `@HD` line, one is added with `VN:1.6`.

This takes time proportional to the size of `header`.

# Examples
```jldoctest
julia> header = XAM2.try_parse_header(ImmutableMemoryView(
           "@HD\\tVN:1.6\\tSO:unsorted\\tGO:query\\n@SQ\\tSN:chr1\\tLN:100\\n"
       ));

julia> header = XAM2.set_sort_order(header, SortOrders.coordinate);

julia> print(String(MemoryView(header)))
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:100
```
"""
function set_sort_order(header::SAMHeader, sort_order::SortOrder)::SAMHeader
    so_field = "SO:" * sort_order_name(sort_order)
    io = IOBuffer()
    print(io, "@HD")
    hd = header.header
    lines = header_lines(header)
    if hd === nothing
        print(io, "\tVN:", LATEST_SAM_VERSION, '\t', so_field)
        pushfirst!(lines, ImmutableMemoryView(take!(io)))
    else
        has_so = false
        for field in split_each(MemoryView(hd), UInt8('\t'))
            # Invariant: Every field is TAG:VALUE, with a two-byte tag
            tag = @inbounds field[1:2]
            if tag == b"SO"
                print(io, '\t', so_field)
                has_so = true
            elseif tag != b"SS" && tag != b"GO"
                write(io, UInt8('\t'), field)
            end
        end
        has_so || print(io, '\t', so_field)
        # Invariant: The @HD line is the first line
        lines[1] = ImmutableMemoryView(take!(io))
    end
    return build_header(lines)
end

"""
    with_references(header::SAMHeader, refs::References)::SAMHeader

Return a copy of `header` where the `@SQ` lines are replaced by lines describing `refs`,
in the same order. Use this to write records with references that differ from those of
the source file, e.g. after subsetting, renaming or reordering references, since writers
require that the `@SQ` lines of the header match the references.

If a reference has the same name and length as an `@SQ` line in `header`, that line is
kept, including its other fields such as `M5` and `AS`. Otherwise, the new line only has
the `SN` and `LN` fields.
The new `@SQ` lines are placed where the first `@SQ` line of `header` was, or after the
`@HD` line if `header` has no `@SQ` lines.

Throw an `ArgumentError` if a reference name is not valid in SAM, see section 1.2.1 of
the SAM specification.

This takes time proportional to the size of `header` plus the size of `refs`.

# Examples
```jldoctest
julia> header = XAM2.try_parse_header(ImmutableMemoryView(
           "@HD\\tVN:1.6\\n@SQ\\tSN:chr1\\tLN:100\\tM5:abc\\n@SQ\\tSN:chr2\\tLN:50\\n"
       ));

julia> header = XAM2.with_references(header, References(["chr3", "chr1"], [10, 100]));

julia> print(String(MemoryView(header)))
@HD	VN:1.6
@SQ	SN:chr3	LN:10
@SQ	SN:chr1	LN:100	M5:abc
```
"""
function with_references(header::SAMHeader, refs::References)::SAMHeader
    old_seqs = Dict{HeaderString, HeaderSQ}(sq.sn => sq for sq in header.seqs)
    sq_lines = Vector{ImmutableMemoryView{UInt8}}(undef, length(refs))
    io = IOBuffer()
    for (i, ref) in zip(eachindex(sq_lines), refs)
        name = ref.name
        if !is_valid_reference_name(codeunits(name))
            throw(ArgumentError(lazy"Reference name \"$name\" is not a valid SAM reference name"))
        end
        print(io, "@SQ\t")
        old = get(old_seqs, name, nothing)
        if old !== nothing && old.ln == ref.length
            write(io, MemoryView(old))
        else
            print(io, "SN:", name, "\tLN:", ref.length)
        end
        sq_lines[i] = ImmutableMemoryView(take!(io))
    end
    lines = header_lines(header)
    sq = str_as_u16("SQ")
    first_sq = findfirst(line -> record_type(line) == sq, lines)
    filter!(line -> record_type(line) != sq, lines)
    if first_sq === nothing
        insert_lines!(lines, sq, sq_lines)
    else
        # The lines before the first @SQ line are unchanged by the filtering
        splice!(lines, first_sq:(first_sq - 1), sq_lines)
    end
    return build_header(lines)
end
