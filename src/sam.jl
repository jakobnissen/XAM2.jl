const SAM_MACHINE = let
    rname_start = re"[0-9A-Za-z!#$%&+./:;?@^_|~-]"
    rname_end = rname_start | re"\*" | re"="
    rname = rname_start * rep(rname_end)

    tab = re"\t"
    unsigned = onall!(re"[0-9]+", :all_number)

    qname = onexit!(re"[!-?A-~]+", :qname)
    flag = onexit!(copy(unsigned), :flag)
    r_name = onenter!(onexit!(re"\*" | rname, :rname), :mark)
    pos = onexit!(copy(unsigned), :pos)
    mapq = onexit!(copy(unsigned), :mapq)
    rnext = onenter!(onexit!(r_name | re"=", :rnext), :mark)
    pnext = onexit!(copy(unsigned), :pnext)
    tlen = onexit!(onenter!(opt("-"), :enter_sign) * copy(unsigned), :tlen)
    cigar = onenter!(re"\*" | re"([0-9]+[MIDNSHP=X])+", :cigar)
    seq = onenter!(re"\*" | re"[A-Za-z=.]+", :seq)
    qual = onenter!(re"[!-~]+", :qual)
    aux = onenter!(re"[ \t!-~]+", :aux)

    line = qname * tab *
        flag * tab *
        r_name * tab *
        pos * tab *
        mapq * tab *
        cigar * tab *
        rnext * tab *
        pnext * tab *
        tlen * tab *
        seq * tab *
        qual *
        opt(tab * aux)

    onexit!(line, :finish)

    Automa.compile(line)
end

const SAM_ACTIONS = Dict{Symbol, Expr}(
    :enter_sign => quote
        is_negative = byte == UInt8('-')
    end,
    :all_number => quote
        (n, overflow) = Base.Checked.mul_with_overflow(n, UInt64(0x0a))
        has_overflowed |= overflow
        (n, overflow) = Base.Checked.add_with_overflow(n, (byte - 0x30) % UInt64)
        has_overflowed |= overflow
    end,
    :mark => quote
        mark = p % UInt32
    end,
    :qname => quote
        # QNAME's length can be 1:254
        p > 255 && return XAMErrors.qname_too_long
        record.qname_len = (p - 1) % UInt8
    end,
    :flag => quote
        n > typemax(UInt16) && return XAMErrors.overflow
        record.flags = new_sam_flags(n % UInt16)
        n = zero(UInt64)
    end,
    :rname => quote
        record.rname_span = (mark, (p - 1) % UInt32)
    end,
    :pos => quote
        n > (typemax(Int32) % UInt64) && return XAMErrors.overflow
        record.pos = n % UInt32
        n = zero(UInt64)
    end,
    :mapq => quote
        n > typemax(UInt8) && return XAMErrors.overflow
        record.mapq = n % UInt8
        n = zero(UInt64)
    end,
    :cigar => quote
        record.cigar_offset = (p - 1) % UInt32
    end,
    :rnext => quote
        record.rnext_span = (mark, (p - 1) % UInt32)
    end,
    :pnext => quote
        n > (typemax(Int32) % UInt64) && return XAMErrors.overflow
        record.pnext = n % UInt32
        n = zero(UInt64)
    end,
    :tlen => quote
        maxval = (typemax(Int32) + is_negative) % UInt64
        n > maxval && return XAMErrors.overflow
        n = n % Int32
        record.tlen = (is_negative ? -n : n)
        n = zero(UInt64)
        is_negative = false
    end,
    :seq => quote
        record.seq_offset = (p - 1) % UInt32
    end,
    :qual => quote
        record.qual_offset = (p - 1) % UInt32
    end,
    :aux => quote
        seen_aux = true
        record.aux_offset = (p - 1) % UInt32
    end,
    :finish => quote
        has_overflowed && return XAMErrors.overflow
    end,
)

mutable struct SAMRecord
    # Data contains the unmodified input line.
    data::Vector{UInt8}
    cigar::Union{Uninit, Nothing, CIGAR, CIGARError}

    rname_span::Tuple{UInt32, UInt32}
    rnext_span::Tuple{UInt32, UInt32}
    cigar_offset::UInt32 # cigar: cigar_offset+1:first(rnext_span)-1
    seq_offset::UInt32 # seq: seq_offset+1:qual_offset-1
    qual_offset::UInt32 # qual: qual_offset+1:aux_offset-1
    aux_offset::UInt32 # aux: aux_offset+1:end
    pos::UInt32
    pnext::UInt32
    tlen::Int32
    flags::SAMFlags
    mapq::UInt8
    qname_len::UInt8

    global function new_uninitialized_sam_record(v::Vector{UInt8})
        return new(v, Uninit())
    end
end

const CTX = Automa.CodeGenContext(; generator = :goto, clean = true)

@eval function parse_record!(record::SAMRecord, data::ImmutableMemoryView{UInt8})
    length(data) > typemax(UInt32) && return XAMErrors.line_too_long
    empty!(record.data)
    append!(record.data, data)
    n = UInt64(0)
    is_negative = false
    has_overflowed = false
    seen_aux = false
    mark = 0
    $(Automa.generate_code(CTX, SAM_MACHINE, SAM_ACTIONS))
    has_overflowed && return XAMErrors.overflow
    if !seen_aux
        record.aux_offset = length(data) % UInt32
    end
    return record
end

SAMRecord() = reset!(new_uninitialized_sam_record(UInt8[]))

const DEFAULT_SAM_DATA_CONTENT = to_memory(collect(b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*"))

const SAM_RECORD_CONTENT = [
    0x05, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, # rname_span: 5:5
    0x0f, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, #  rnext_span: 15:15
    0x0c, 0x00, 0x00, 0x00, # cigar offset: 12
    0x14, 0x00, 0x00, 0x00, # seq offset: 20
    0x16, 0x00, 0x00, 0x00, # qual offset: 22
    0x17, 0x00, 0x00, 0x00, # aux offset: 23 (when missing)
    0x00, 0x00, 0x00, 0x00, # pos
    0x00, 0x00, 0x00, 0x00, # pnext
    0x00, 0x00, 0x00, 0x00, # tlen
    0x04, 0x00, # flags: unmapped
    0xff, # mapq: unavailable,
    0x00, # qname_len
] |> to_memory


const FIXED_SAM_OFFSET = fieldoffset(SAMRecord, 3)
const FIXED_SAM_SIZE = 48

@assert sizeof(SAMRecord) == FIXED_SAM_SIZE + FIXED_SAM_OFFSET
@assert sizeof(SAM_RECORD_CONTENT) == FIXED_SAM_SIZE

function reset_fixed_fields!(record::SAMRecord)
    GC.@preserve record SAM_RECORD_CONTENT begin
        dst = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_SAM_OFFSET
        src = pointer(SAM_RECORD_CONTENT)
        @ccall memcpy(dst::Ptr{UInt8}, src::Ptr{UInt8}, Csize_t(FIXED_SAM_SIZE)::Csize_t)::Cvoid
    end
    return record
end

function reset!(record::SAMRecord)
    resize!(record.data, length(DEFAULT_SAM_DATA_CONTENT))
    @inbounds copyto!(MemoryView(record.data), ImmutableMemoryView(DEFAULT_SAM_DATA_CONTENT))
    record.cigar = Uninit()
    return reset_fixed_fields!(record)
end

function Base.copy(record::SAMRecord)
    vector_copy = copy(record.data)
    new = new_uninitialized_sam_record(vector_copy)
    GC.@preserve record new begin
        dst = Ptr{UInt8}(pointer_from_objref(new)) + FIXED_SAM_OFFSET
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_SAM_OFFSET
        @ccall memcpy(dst::Ptr{UInt8}, src::Ptr{UInt8}, Csize_t(FIXED_SAM_SIZE)::Csize_t)::Cvoid
    end
    record.cigar = copy_cigar(record.cigar, record.data, vector_copy)
    return new
end


# Short form: Mapped, read name, strand
function Base.show(io::IO, record::SAMRecord)
    mapped = is_mapped(record) ? "Mapped" : "Unmapped"
    segment = is_forward(record) ? "Fw" : (is_reverse(record) ? "Rv" : "Md")
    qname = query_name(record)
    qn = if qname === nothing
        ""
    else
        ", " * repr(qname)
    end
    return print(io, "SAMRecord<", mapped, ", ", segment, qn, '>')
end

function Base.show(io::IO, ::MIME"text/plain", record::SAMRecord)
    rtruncate = Base.rtruncate

    buf = IOBuffer()
    print(buf, is_mapped(record) ? "Mapped " : "Unmapped ")
    summary(buf, record)
    width = max(displaysize()[2] - 8, 20)
    println(buf, ':')

    # Qname
    print(buf, "  Name: ")
    qname = query_name(record)
    if isnothing(qname)
        println(buf, "< missing >")
    else
        println(buf, rtruncate(repr(qname), width - 12))
    end

    # Flags
    flagstr = join([string(i) for i in flags(record)], ", ")
    println(buf, "  Flags: ", repr(Integer(flags(record))), ": ", Base.rtruncate(flagstr, width - 12))

    # Ref name
    print(buf, "  Ref: ")
    ref = reference_name(record)
    if ref === nothing
        println(buf, "< missing >")
    else
        println(buf, rtruncate(ref, width - 12))
    end

    # Mapping span
    print(buf, "  Ref pos: ")
    span = mapping_span(record)
    if span !== nothing
        println(buf, span)
    else
        left = left_position(record)
        if left !== nothing
            println(buf, left)
        else
            println(buf, "< missing >")
        end
    end

    # Seq
    print(buf, "  Seq: ")
    seq = sequence(record)
    if seq === nothing
        println(buf, "< missing >")
    else
        println(buf, rtruncate(string(seq), width - 12))
    end

    # CIGAR
    print(buf, "  CIGAR: ")
    cig = cigar(record)
    if cig === nothing
        println(buf, "< missing >")
    elseif cig isa CIGARError
        println(buf, "< malformed CIGAR >")
    else
        str = String(MemoryView(CIGAR(cig)))
        println(buf, rtruncate(str, width - 12))
    end

    # Auxiliary data
    aux = aux_data(record)
    if !isempty(aux)
        println(buf, "  Auxiliary:")
        auxbuf = IOBuffer()
        show(auxbuf, MIME"text/plain"(), aux_data(record))
        seekstart(auxbuf)
        readline(auxbuf)
        for line in eachline(auxbuf)
            println(buf, "    ", line)
        end
    end
    v = take!(buf)
    pop!(v) # trailing newline
    return write(io, v)
end


function query_name(record::SAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    view = view[1:record.qname_len]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = ImmutableMemoryView(record.data)
    span = record.rname_span
    view = view[first(span):last(span)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function cigar(record::SAMRecord)::Union{Nothing, CIGAR, CIGARError}
    record.cigar === Uninit() || return record.cigar
    view = ImmutableMemoryView(record.data)
    view = view[(record.cigar_offset + 1):(first(record.rname_span) - 1)]
    length(view) == 1 && @inbounds view[1] == UInt8('*') && return nothing
    return record.cigar = CIGARStrings.try_parse(CIGAR, view)
end

function next_reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = ImmutableMemoryView(record.data)
    span = record.rnext_span
    view = view[first(span):last(span)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function sequence(record::SAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    view = view[(record.seq_offset + 1):(record.qual_offset - 1)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

# TODO: This returns a vector with offset + 33 whereas the BAM one does not.
# should they be the same?
function quality(record::SAMRecord)::Union{Nothing, AbstractString}
    view = ImmutableMemoryView(record.data)
    view = view[(record.qual_offset + 1):(record.aux_offset - 1)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function aux_data(record::SAMRecord)
    return XAMAuxData.SAM.Auxiliary(record.data, record.aux_offset + 1)
end

function left_position(record::SAMRecord)::Union{Int, Nothing}
    result = record.pos
    return iszero(result) ? nothing : result
end

mutable struct SAMReader{T <: AbstractBufReader}
    # TODO: Must be able to buffer a line - put in docstring
    const reader::T
    const record::SAMRecord
    const offset_first_record::Int
    const copy::Bool
    error::Union{Nothing, XAMError}
end

function SAMReader(reader::AbstractBufReader; copy::Bool = true)
    (header, offset_first_record) = read_initial_sam_header!(reader)
    record = SAMRecord()
    return SAMReader{typeof(reader)}(reader, record, offset_first_record, header, copy, nothing)
end

function SAMReader(io::IO; copy::Bool = true)
    return SAMReader(BufReader(io); copy)
end

function Base.show(io::IO, reader::SAMReader)
    print(io, '<')
    summary(io, reader)
    return print(io, '>')
end

function Base.show(io::IO, ::MIME"text/plain", reader::SAMReader)
    summary(io, reader)
    print(io, ':')
    return print(io, "TODO!")
end


function read_initial_sam_header!(reader::AbstractBufReader)
    # We need to store a copy of the underlying data, and can't just store
    # a view into the reader, because we mutate the reader
    # TODO: Is it too inefficient to do two passes over the header, one to
    # find the end of the buffer, and one to parse?
    # Can we do this in one pass?
    buffer = VecWriter(undef, 4096)
    for line_view in line_views(reader; chomp = false)
        if isempty(line_view) || @inbounds(line_view[1]) != UInt8('@')
            break
        end
        write(buffer, line_view)
    end
    data = get_unflushed(buffer)
    offset = length(data)
    return (parse_sam_header(ImmutableMemoryView(data)), offset)
end

function seek_first_record(reader::SAMReader)
    seek(reader.reader, reader.offset_first_record)
    reset!(reader.record)
    return reader
end

function try_read!(reader::SAMReader)::Union{SAMRecord, Nothing, XAMError}
    # Get a view of the current line, stripped from newline
    isnothing(reader.error) || return reader.error
    line_view = get_next_line_view(reader.reader)
    line_view === nothing && return nothing
    stripped = _chomp(line_view)
    result = parse_record!(reader.record, stripped)
    @inbounds consume(reader.reader, length(line_view))
    return result
end

function get_next_line_view(
        reader::AbstractBufReader
    )::Union{ImmutableMemoryView{UInt8}, Nothing}
    scan_from = 1
    buffer = get_buffer(reader)::ImmutableMemoryView{UInt8}
    while true
        pos = findnext(==(0x0a), buffer, scan_from)
        if pos !== nothing
            return @inbounds(buffer[1:pos])
        end
        scan_from = length(buffer) + 1
        n_filled = fill_buffer(reader)
        if n_filled === nothing
            error(
                "Could not store full line in buffer. " *
                    "SAMReader must be constructed from an `AbstractBufReader` that can store " *
                    "an entire line in its buffer."
            )
        end
        new_buffer = get_buffer(reader)::ImmutableMemoryView{UInt8}
        if length(new_buffer) ≤ length(buffer)
            return isempty(new_buffer) ? nothing : new_buffer
        end
        buffer = new_buffer
    end
    return # unreachable
end

function _chomp(x::ImmutableMemoryView{UInt8})::ImmutableMemoryView{UInt8}
    len = if isempty(x)
        0
    else
        has_lf = @inbounds(x[end]) == 0x0a
        two_bytes = length(x) > 1
        has_cr = has_lf & two_bytes & (@inbounds(x[length(x) - two_bytes]) == 0x0d)
        length(x) - (has_lf + has_cr)
    end
    @inbounds return x[1:len]
end
