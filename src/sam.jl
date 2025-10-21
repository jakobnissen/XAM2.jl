# TODO: Instead of packing fields, store the entire line directly?
mutable struct SAMRecord
    # QNAME, RNAME, CIGAR, RNEXT, SEQ, QUAL, AUX, in order
    data::Vector{UInt8}
    cigar::Union{Uninit, Nothing, CIGAR, CIGARError}
    rname_offset::UInt32
    cigar_offset::UInt32
    rnext_offset::UInt32
    seq_offset::UInt32
    qual_offset::UInt32
    aux_offset::UInt32
    pos::Int32
    pnext::Int32
    tlen::Int32
    flags::SAMFlags
    mapq::UInt8
    _padding::UInt8

    global function new_uninitialized_sam_record(v::Vector{UInt8})
        return new(v, Uninit())
    end
end

SAMRecord() = reset!(new_uninitialized_sam_record(UInt8[]))

const SAM_RECORD_CONTENT = [
    0x01, 0x00, 0x00, 0x00, # rname offset: query name is '*'
    0x02, 0x00, 0x00, 0x00, # cigar offset: rname is '*'
    0x03, 0x00, 0x00, 0x00, # rnext offset: cigar is '*'
    0x04, 0x00, 0x00, 0x00, # seq offset: rnext is '*'
    0x05, 0x00, 0x00, 0x00, # qual offset: seq is '*'
    0x06, 0x00, 0x00, 0x00, # aux offset: qual is '*'
    0x00, 0x00, 0x00, 0x00, # pos
    0x00, 0x00, 0x00, 0x00, # pnext
    0x00, 0x00, 0x00, 0x00, # tlen
    0x04, 0x00, # flags: unmapped
    0xff, # mapq: unavailable,
    0x00, # padding: zeroed out
] |> to_memory

const FIXED_SAM_OFFSET = fieldoffset(SAMRecord, 3)
const FIXED_SAM_SIZE = 40

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
    resize!(record.data, 6)
    copyto!(record.data, b"******")
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
    view = view[1:record.rname_offset]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = ImmutableMemoryView(record.data)
    view = view[(record.rname_offset + 1):record.cigar_offset]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function cigar(record::SAMRecord)::Union{Nothing, CIGAR, CIGARError}
    record.cigar === Uninit() || return record.cigar
    view = ImmutableMemoryView(record.data)
    view = view[record.cigar_offset:record.rnext_offset]
    length(view) == 1 && @inbounds view[1] == UInt8('*') && return nothing
    return record.cigar = CIGARStrings.try_parse(CIGAR, view)
end

function next_reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = ImmutableMemoryView(record.data)
    view = view[(record.rnext_offset + 1):record.seq_offset]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function sequence(record::SAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    view = view[(record.seq_offset + 1):record.qual_offset]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

# TODO: This returns a vector with offset + 33 whereas the BAM one does not.
function quality(record::SAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    view = ImmutableMemoryView(record.data)
    view = view[(record.qual_offset + 1):record.aux_offset]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        view
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
    const header::SAMHeader
    const fields::Memory{UnitRange{UInt32}}
    const copy::Bool
    error::Union{Nothing, XAMError}
end

function SAMReader(reader::AbstractBufReader; copy::Bool = true)
    (header, offset_first_record) = read_initial_sam_header!(reader)
    record = SAMRecord()
    fields = Memory{UnitRange{UInt32}}(undef, 12)
    @inbounds for i in eachindex(fields)
        fields[i] = UInt32(1):UInt32(0)
    end
    return SAMReader{typeof(reader)}(reader, record, offset_first_record, header, fields, copy, nothing)
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

# TODO: Should we verify the content of the fields, e.g. that qname
# only contains the right bytes?
# TODO: Lots of boundschecks can be elided, review after testing
function try_read!(reader::SAMReader)::Union{SAMRecord, Nothing, XAMError}
    # Get a view of the current line, stripped from newline
    isnothing(reader.error) || return reader.error
    line_view = get_next_line_view(reader.reader)
    line_view === nothing && return nothing
    stripped = _chomp(line_view)

    # Read each tab-separated field
    fields = reader.fields
    maybe_error = load_fields(fields, stripped)
    mem = ImmutableMemoryView(parent(stripped))
    maybe_error isa XAMError && return maybe_error

    # Ensure enough room in record's data, and get a view of the
    # data to write into
    record = reader.record
    data = record.data
    length(data) < length(line_view) && resize!(data, length(line_view))
    remaining_data = MemoryView(data)
    n_written = 0

    # qname
    field = @inbounds(mem[fields[1]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.rname_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # flags
    flagint = tryparse_base10(UInt16, @inbounds(mem[fields[2]]))
    flagint === nothing && return XAMErrors.invalid_integer
    # constructors ignores top 4 bits per specs
    record.flags = SAMFlags(flagint)

    # rname
    field = @inbounds(mem[fields[3]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.cigar_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # pos
    posint = tryparse_base10(UInt32, @inbounds(mem[fields[4]]))
    posint === nothing && return XAMErrors.invalid_integer
    posint > (typemax(Int32) % UInt32) && return XAMErrors.invalid_integer
    record.pos = posint % Int32

    # mapq
    mapqint = tryparse_base10(UInt8, @inbounds(mem[fields[5]]))
    mapqint === nothing && return XAMErrors.invalid_integer
    record.mapq = mapqint

    # cigar
    field = @inbounds(mem[fields[6]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.rnext_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # rnext
    field = @inbounds(mem[fields[7]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.seq_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # pnext
    pnextint = tryparse_base10(UInt32, @inbounds(mem[fields[8]]))
    pnextint === nothing && return XAMErrors.invalid_integer
    pnextint > (typemax(Int32) % UInt32) && return XAMErrors.invalid_integer
    record.pnext = pnextint % Int32

    # tlen
    pnextint = tryparse_base10(Int32, @inbounds(mem[fields[9]]))
    pnextint === nothing && return XAMErrors.invalid_integer
    record.tlen = pnextint

    # seq
    field = @inbounds(mem[fields[10]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.qual_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # qual
    field = @inbounds(mem[fields[11]])
    @inbounds copyto!(remaining_data, field)
    n_written += length(field)
    record.aux_offset = n_written % UInt32
    remaining_data = @inbounds remaining_data[(length(field) + 1):end]

    # Aux, if present
    aux_span = @inbounds fields[12]
    if !isempty(aux_span)
        aux = @inbounds mem[aux_span]
        @inbounds copyto!(remaining_data, aux)
        n_written += length(aux)
    end

    resize!(record.data, n_written)
    @inbounds consume(reader.reader, length(line_view))

    return record
end

function load_fields(
        mem::Memory{UnitRange{UInt32}},
        stripped::ImmutableMemoryView{UInt8},
    )::Union{Nothing, XAMError}
    length(mem) == 12 || error() # cheaper assert
    i = 0
    for field in split_each(stripped, UInt8('\t'))
        i += 1
        span = only(parentindices(field))
        @inbounds mem[i] = (first(span) % UInt32):(last(span) % UInt32)
        i == 11 && break
    end
    i < 11 && return XAMErrors.missing_fields
    qual_span = @inbounds mem[11]
    aux_start = last(qual_span) + UInt32(2)
    aux_stop = last(only(parentindices(stripped))) % UInt32
    @inbounds mem[12] = aux_start:aux_stop
    return nothing
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
