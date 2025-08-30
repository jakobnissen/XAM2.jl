struct HeaderTag
    a::UInt8
    b::UInt8

    global new_tag(a::UInt8, b::UInt8) = new(a, b)
end

try_tag(x::AbstractString) = try_tag(codeunits(x))

function try_tag(x::AbstractVector{UInt8})
    length(x) == 2 || return nothing
    a = @inbounds x[1]
    b = @inbounds x[2]
    upper = UInt8('A'):UInt8('Z')
    lower = UInt8('a'):UInt8('z')
    digit = UInt8('0'):UInt8('9')
    if (
            (!in(a, upper) & !in(a, lower)) |
                (!in(b, upper) & !in(b, lower) & !in(b, digit))
        )
        return nothing
    end
    return new_tag(a, b)
end

const TAG_HD = try_tag("HD")
const TAG_SQ = try_tag("SQ")
const TAG_RG = try_tag("RG")
const TAG_PG = try_tag("PG")
const TAG_CO = try_tag("CO")


struct HeaderField
    tag::HeaderTag
    data::ImmutableMemoryView{UInt8}
end

struct MetaInfo
    tag::HeaderTag
    fields::Vector{HeaderField}
end

struct SAMHeader
    # Only one header, and it must be first if present
    metainfo::Vector{MetaInfo}
end

function parse_sam_header(data::ImmutableMemoryView{UInt8})
    v = MetaInfo[]
    for (line_index, line) in enumerate(eachline(CursorReader(data)))
        mi = parse_metainfo(line)
        mi.tag == TAG_HD && line_index > 1 && error("TODO!")
        push!(v, mi)
    end
    return SAMHeader(v)
end

# TODO: Do not error, return error value
function parse_metainfo(data::ImmutableMemoryView{UInt8})::MetaInfo
    # Assumed must not contain a newline
    length(data) < 3 && error("Truncated")
    @inbounds data[1] == UInt8('@') || error("Not @")

    tag = try_tag(@inbounds data[2], @inbounds data[3])
    tag === nothing && error("Invalid header tag")

    metainfo = MetaInfo(tag, HeaderField[])

    rest = if length(data) == 3
        @inbounds data[3:end]
    else
        @inbounds data[4] == UInt8('\t') || error("Bad tag")
        @inbounds data[5:end]
    end

    # TODO: This is private - make public, or maybe move to MemoryViews?
    fields = XAMAuxData.DelimitedIterator(rest, UInt8('\t'))

    return if tag == TAG_HD
        parse_header_hd(metainfo, fields)
    elseif tag == TAG_SQ
        parse_header_sq(metainfo, fields)
    elseif tag == TAG_RG
        parse_header_rg(metainfo, fields)
    elseif tag == TAG_PG
        parse_header_pg(metainfo, fields)
    elseif tag == TAG_CO
        parse_header_co(metainfo, fields)
    else
        error("Invalid header tag")
    end
end
