const SAM_MACHINE = let
    # See the SAM specification, section 1.4
    tab = re"\t"
    digits = onall!(re"[0-9]+", :digit)
    reference_name = re"\*|[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*"

    qname = onexit!(re"[!-?A-~]+", :qname)
    flag = onexit!(copy(digits), :flag)
    rname = onenter!(onexit!(copy(reference_name), :rname), :mark)
    pos = onexit!(copy(digits), :pos)
    mapq = onexit!(copy(digits), :mapq)
    cigar = onenter!(re"\*|([0-9]+[MIDNSHP=X])+", :cigar)
    rnext = onenter!(onexit!(copy(reference_name) | re"=", :rnext), :mark)
    pnext = onexit!(copy(digits), :pnext)
    tlen = onexit!(opt(onenter!(re"-", :negative) | re"\+") * copy(digits), :tlen)
    seq = onenter!(re"\*|[A-Za-z=.]+", :seq)
    qual = onenter!(re"[!-~]+", :qual)
    # Auxiliary data is validated lazily
    aux = onenter!(re"[ \t!-~]+", :aux)

    line = qname * tab *
        flag * tab *
        rname * tab *
        pos * tab *
        mapq * tab *
        cigar * tab *
        rnext * tab *
        pnext * tab *
        tlen * tab *
        seq * tab *
        qual *
        opt(tab * aux)

    Automa.compile(line)
end

# Numeric fields are accumulated in `n`, and checked when the field is exited.
# Out-of-range numbers return the error of the field. When a field is exited, `p` is
# the index of the byte after it.
const SAM_ACTIONS = Dict{Symbol, Expr}(
    :digit => quote
        (n, overflow) = Base.Checked.mul_with_overflow(n, UInt64(0x0a))
        has_overflowed |= overflow
        (n, overflow) = Base.Checked.add_with_overflow(n, (byte - 0x30) % UInt64)
        has_overflowed |= overflow
    end,
    :mark => :(mark = p),
    :qname => quote
        # QNAME's length can be 1:254
        p > 255 && return sam_field_error(data, XAMErrors.invalid_read_name, p)
        record.qname_len = (p - 1) % UInt8
    end,
    :flag => quote
        (has_overflowed | (n > typemax(UInt16))) && return sam_field_error(data, XAMErrors.invalid_flags, p)
        # Clears undefined flag bits
        record.flags = SAMFlags(n % UInt16)
        n = zero(UInt64)
    end,
    :rname => :(record.rname_span = (mark % UInt32, (p - 1) % UInt32)),
    :pos => quote
        (has_overflowed | (n > (typemax(Int32) % UInt64))) && return sam_field_error(data, XAMErrors.invalid_position, p)
        record.pos = n % UInt32
        n = zero(UInt64)
    end,
    :mapq => quote
        (has_overflowed | (n > typemax(UInt8))) && return sam_field_error(data, XAMErrors.invalid_mapq, p)
        record.mapq = n % UInt8
        n = zero(UInt64)
    end,
    :cigar => :(record.cigar_offset = (p - 1) % UInt32),
    :rnext => :(record.rnext_span = (mark % UInt32, (p - 1) % UInt32)),
    :pnext => quote
        (has_overflowed | (n > (typemax(Int32) % UInt64))) &&
            return sam_field_error(data, XAMErrors.invalid_next_position, p)
        record.pnext = n % UInt32
        n = zero(UInt64)
    end,
    :negative => :(is_negative = true),
    :tlen => quote
        maxval = (typemax(Int32) % UInt64) + is_negative
        (has_overflowed | (n > maxval)) && return sam_field_error(data, XAMErrors.invalid_template_length, p)
        n32 = n % Int32
        record.tlen = is_negative ? -n32 : n32
    end,
    :seq => :(record.seq_offset = (p - 1) % UInt32),
    :qual => :(record.qual_offset = (p - 1) % UInt32),
    :aux => quote
        seen_aux = true
        record.aux_offset = (p - 1) % UInt32
    end,
)

"""
    SAMRecord

A `SAMRecord` represents a non-header line of a SAM file.

This type is semi-lazily validated; some fields, such as the sequence, quality and CIGAR
will not be validated until queries.

Construct a valid, but empty record with `SAMRecord()`.

!!! warning
    Mutating its backing memory will invalidate the record and
    lead to logic bugs and garbage results, but should never cause
    undefined behaviour in the absence of data races

# Examples
```jldoctest
julia> empty_rec = SAMRecord();

julia> is_mapped(empty_rec)
false

julia> query_name(empty_rec) === nothing
true

julia> aux_data(empty_rec) |> isempty
true
```
"""
mutable struct SAMRecord
    # Structural invariants. Every `SAMRecord` satisfies these, so accessors may rely on them,
    # e.g. using `@inbounds`. Records from `try_parse!` are checked by `SAM_MACHINE`, and
    # records created any other way are valid by construction.
    # The fields of `SAMRecord` are internal, so users can't violate these.
    # 1. `length(data) < typemax(UInt32)`
    # 2. `data` is a SAM line without the trailing newline, with 11 or more tab-separated
    #    fields, and the offsets below point to the correct fields, as described in
    #    their comments.
    # 3. If SEQ is `*`, QUAL is `*`, and otherwise, QUAL is `*` or has the same length as SEQ.
    # Records from `try_parse!` additionally have syntactically valid fields, except the
    # auxiliary data, which is validated when accessed, and the numeric fields are stored
    # in this struct. Records from `copy!(::SAMRecord, ::BAMRecord)` may have invalid
    # content copied from the BAM record, e.g. out-of-range qualities.
    # The CIGAR's content (e.g. that it matches the sequence length) is validated
    # when accessed.

    # Data contains the unmodified input line.
    data::Vector{UInt8}
    cigar::Union{Uninit, Nothing, CIGAR, CIGARError, XAMError}

    # Inclusive start and stop indices of RNAME and RNEXT in `data`
    rname_span::Tuple{UInt32, UInt32}
    rnext_span::Tuple{UInt32, UInt32}
    cigar_offset::UInt32 # cigar: cigar_offset+1:first(rnext_span)-2
    seq_offset::UInt32 # seq: seq_offset+1:qual_offset-1
    qual_offset::UInt32 # qual: qual_offset+1:aux_offset-1
    # Index of the tab before the auxiliary data, or `length(data) + 1` if there
    # is no auxiliary data. Aux: aux_offset+1:end
    aux_offset::UInt32
    pos::UInt32 # 1-based, 0 if missing
    pnext::UInt32 # 1-based, 0 if missing
    tlen::Int32
    flags::SAMFlags
    mapq::UInt8
    qname_len::UInt8

    global function new_uninitialized_sam_record(v::Vector{UInt8})
        return new(v, Uninit())
    end
end

"""
    XAM2.try_parse!(
        record::SAMRecord,
        data::ImmutableMemoryView{UInt8}
    )::Union{SAMRecord, XAMError}

Parse the line `data` (without a trailing newline) as a SAM record into `record`,
returning it. The structure of the record is validated, e.g. that it has 11 fields,
and that the fixed-length fields are valid. The CIGAR and auxiliary data are
validated when accessed.

If the data is invalid, `record` is reset to an empty record, and an `XAMError` is
returned. Its kind identifies which field is invalid, and its column is the position
of the invalid byte, or the start of the field if its value is out of range.

This function does not throw if `data` is malformed.
"""
function try_parse!(record::SAMRecord, data::ImmutableMemoryView{UInt8})::Union{SAMRecord, XAMError}
    err = if length(data) ≥ typemax(UInt32)
        XAMError(XAMErrors.line_too_long)
    else
        copyto!(resize!(record.data, length(data)), data)
        record.cigar = Uninit()
        parse_sam_fields!(record, data)
    end
    if err === nothing
        err = check_seq_qual(record)
    end
    if err !== nothing
        # The record's fields may be partially overwritten, breaking its invariants
        reset!(record)
        return err
    end
    return record
end

@eval function parse_sam_fields!(record::SAMRecord, data::ImmutableMemoryView{UInt8})::Union{Nothing, XAMError}
    n = UInt64(0)
    is_negative = false
    has_overflowed = false
    seen_aux = false
    mark = 0
    $(Automa.generate_init_code(CTX, SAM_MACHINE))
    $(Automa.generate_exec_code(CTX, SAM_MACHINE, SAM_ACTIONS))
    # cs is 0 on success, negative on an invalid byte at `p`, and positive if
    # the input ended before the line was complete.
    iszero(cs) || return sam_syntax_error(data, p, cs)
    if !seen_aux
        record.aux_offset = (length(data) + 1) % UInt32
    end
    return nothing
end

# The error of each of the 11 mandatory fields, in order
const SAM_FIELD_ERRORS = (
    XAMErrors.invalid_read_name,
    XAMErrors.invalid_flags,
    XAMErrors.invalid_reference_name,
    XAMErrors.invalid_position,
    XAMErrors.invalid_mapq,
    XAMErrors.invalid_cigar,
    XAMErrors.invalid_next_reference_name,
    XAMErrors.invalid_next_position,
    XAMErrors.invalid_template_length,
    XAMErrors.invalid_sequence,
    XAMErrors.invalid_quality,
)

# Slow path, only run when the machine fails: Find the field where the error is.
@noinline function sam_syntax_error(data::ImmutableMemoryView{UInt8}, p::Int, cs::Int)::XAMError
    # On error, p is the index of the invalid byte, and if the input ended early,
    # p is length(data) + 1.
    field = count(==(UInt8('\t')), @inbounds(data[1:(p - 1)])) + 1
    kind = if cs > 0 && field < 11
        # Input ended before the 11th field
        XAMErrors.missing_fields
    elseif field > 11
        XAMErrors.invalid_aux_data
    else
        @inbounds(SAM_FIELD_ERRORS[field])
    end
    return XAMError(kind; column = p)
end

# Slow path for a field with an out-of-range value, where `p` is the index after the
# field. The column of the error is the start of the field.
@noinline function sam_field_error(data::ImmutableMemoryView{UInt8}, kind::XAMErrorKind, p::Int)::XAMError
    tab = findprev(==(UInt8('\t')), data, p - 1)
    return XAMError(kind; column = something(tab, 0) + 1)
end

# Check that QUAL is missing if SEQ is, and else has the same length.
function check_seq_qual(record::SAMRecord)::Union{Nothing, XAMError}
    seq = seq_view(record)
    qual = qual_view(record)
    is_missing(qual) && return nothing
    kind = if is_missing(seq)
        XAMErrors.qual_no_seq
    elseif length(seq) != length(qual)
        XAMErrors.qual_seq_mismatch
    else
        return nothing
    end
    return XAMError(kind; column = (record.qual_offset % Int) + 1)
end

# Check if a field is `*`, meaning missing
is_missing(mem::ImmutableMemoryView{UInt8}) = length(mem) == 1 && @inbounds(mem[1]) == UInt8('*')

# Invariant: The offsets point to the correct fields, see the definition of `SAMRecord`.
# Computed as `Int` to avoid overflow.
function seq_view(record::SAMRecord)
    span = ((record.seq_offset % Int) + 1):((record.qual_offset % Int) - 1)
    return @inbounds ImmutableMemoryView(record.data)[span]
end

function qual_view(record::SAMRecord)
    span = ((record.qual_offset % Int) + 1):((record.aux_offset % Int) - 1)
    return @inbounds ImmutableMemoryView(record.data)[span]
end

function cigar_view(record::SAMRecord)
    span = ((record.cigar_offset % Int) + 1):((first(record.rnext_span) % Int) - 2)
    return @inbounds ImmutableMemoryView(record.data)[span]
end

function span_view(record::SAMRecord, span::Tuple{UInt32, UInt32})
    return @inbounds ImmutableMemoryView(record.data)[(first(span) % Int):(last(span) % Int)]
end

"""
    parse!(record::SAMRecord, data::ImmutableMemoryView{UInt8})::SAMRecord

Parse the `data` into a `record`, or throw a `XAMError` if the data is invalid.
"""
function parse!(record::SAMRecord, data::ImmutableMemoryView{UInt8})
    result = try_parse!(record, data)
    result isa XAMError && throw(result)
    return result
end

SAMRecord() = reset!(new_uninitialized_sam_record(UInt8[]))

const DEFAULT_SAM_DATA_CONTENT = Tuple(b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*")

const SAM_RECORD_CONTENT = [
    0x05, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, # rname_span: 5:5
    0x0f, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, #  rnext_span: 15:15
    0x0c, 0x00, 0x00, 0x00, # cigar offset: 12
    0x14, 0x00, 0x00, 0x00, # seq offset: 20
    0x16, 0x00, 0x00, 0x00, # qual offset: 22
    0x18, 0x00, 0x00, 0x00, # aux offset: 24 (when missing, length of data + 1)
    0x00, 0x00, 0x00, 0x00, # pos
    0x00, 0x00, 0x00, 0x00, # pnext
    0x00, 0x00, 0x00, 0x00, # tlen
    0x04, 0x00, # flags: unmapped
    0xff, # mapq: unavailable,
    0x01, # qname_len
] |> Tuple


const FIXED_SAM_OFFSET = fieldoffset(SAMRecord, 3)
const FIXED_SAM_SIZE = 48

@assert sizeof(SAMRecord) == FIXED_SAM_SIZE + FIXED_SAM_OFFSET
@assert sizeof(SAM_RECORD_CONTENT) == FIXED_SAM_SIZE

function reset_fixed_fields!(record::SAMRecord)
    # See comment in `reset_fixed_fields!(::BAMRecord)` for why this is a tuple
    GC.@preserve record begin
        dst = Ptr{typeof(SAM_RECORD_CONTENT)}(pointer_from_objref(record) + FIXED_SAM_OFFSET)
        unsafe_store!(dst, SAM_RECORD_CONTENT)
    end
    return record
end

"""
    reset!(::SAMRecord)

Resets the `SAMRecord` to its default state, deleting all information from it.
After this function, the record is equal to a `SAMRecord()`.
"""
function reset!(record::SAMRecord)
    data = record.data
    resize!(data, length(DEFAULT_SAM_DATA_CONTENT))
    for i in 1:length(DEFAULT_SAM_DATA_CONTENT)
        @inbounds data[i] = @inbounds(DEFAULT_SAM_DATA_CONTENT[i])
    end
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
    new.cigar = copy_cigar!(record.cigar, record.data, vector_copy)
    return new
end


# Short form: Mapped, read name, strand
function Base.show(io::IO, record::SAMRecord)
    mapped = is_mapped(record) ? "Mapped" : "Unmapped"
    segment = is_first(record) ? "First" : (is_last(record) ? "Last" : "Middle")
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
    if span isa UnitRange{Int}
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
    elseif cig isa Union{CIGARError, XAMError}
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
    # Invariant: QNAME is the first field, and has length `qname_len`
    view = @inbounds ImmutableMemoryView(record.data)[1:record.qname_len]
    return is_missing(view) ? nothing : StringView(view)
end

function reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = span_view(record, record.rname_span)
    return is_missing(view) ? nothing : StringView(view)
end

function cigar(record::SAMRecord)::Union{Nothing, CIGAR, CIGARError, XAMError}
    record.cigar === Uninit() || return record.cigar
    view = cigar_view(record)
    if is_missing(view)
        record.cigar = nothing
        return nothing
    end
    cig = CIGARStrings.try_parse(CIGAR, view)
    if !isa(cig, CIGARError)
        seq = sequence(record)
        seq_len = seq === nothing ? 0 : ncodeunits(seq)
        cig = check_query_length(cig, seq_len)
    end
    return record.cigar = cig
end

function next_reference_name(record::SAMRecord)::Union{Nothing, AbstractString}
    view = span_view(record, record.rnext_span)
    return if is_missing(view)
        nothing
    elseif length(view) == 1 && @inbounds view[1] == UInt8('=')
        reference_name(record)
    else
        StringView(view)
    end
end

"""
    sequence(record::SAMRecord)::Union{AbstractString, Nothing}

Get the query sequence of `rec` as an `AbstractString`, or
`nothing` if the sequence is not available.

The sequence is a UTF-8 encoded view into the record.

!!! warning
    The returned sequence aliases the record.
    Mutating the record will invalidate the sequence.
"""
function sequence(record::SAMRecord)::Union{AbstractString, Nothing}
    view = seq_view(record)
    return is_missing(view) ? nothing : StringView(view)
end

"""
    quality_string(record::SAMRecord)::Union{Nothing, AbstractString}

Return the quality as a PHRED+33-encoded `AbstractString`,
or `nothing` if the quality is not available.

!!! warning
    The returned string aliases the record.
    Mutating the record will invalidate the quality string.
"""
function quality_string(record::SAMRecord)::Union{Nothing, AbstractString}
    view = qual_view(record)
    return is_missing(view) ? nothing : StringView(view)
end

"""
    aux_data(record::SAMRecord)::XAMAuxData.SAM.MutableAuxiliary

Return the auxiliary data as a mutable `XAMAuxData.SAM.Auxiliary`.
If the record does not contain any auxiliary fields, the returned object
is empty.

!!! warning
    The returned `Auxiliary` object aliases the record.
    Mutating the record will invalidate the auxiiary data.

See the documentation of `XAMAuxData.SAM.Auxiliary` for more details.
"""
function aux_data(record::SAMRecord)
    # If there is no auxiliary data, `aux_offset + 1` is `length(data) + 2`
    start = min((record.aux_offset % Int) + 1, length(record.data) + 1)
    return XAMAuxData.SAM.Auxiliary(record.data, start)
end

function left_position(record::SAMRecord)::Union{Int, Nothing}
    result = record.pos % Int
    return iszero(result) ? nothing : result
end

"""
    mutable struct SAMReader{T <: AbstractBufReader}
    
Lazy reader of SAM files.

`rdr = SAMReader` is an iterator over `SAMRecord`.
The header can be obtained with [`get_header`](@ref)

    SAMReader(
        io::Union{IO, AbstractBufReader};
        copy::Bool = true,
    )::SAMReader

Construct a `SAMReader` from `io`, reading and validating the header.
This reader takes control of `io`. Mutating `io` after the reader is constructed
is not permitted. Closing the resulting reader will close `io`.
If the header is malformed, `io` is closed and an `XAMError` is thrown.
See [`XAM2.try_init`](@ref) for a version that returns the error instead.

If `copy` is `false`, every record returned by iterating the reader or from
`XAM2.try_read!` will be the same record, overwritten in-place. This improves
efficiency, but new records will invalidate every previous one.
To prevent aliasing, `copy` can be called manually on records.
"""
mutable struct SAMReader{T <: AbstractBufReader}
    # TODO: Must be able to buffer a line - put in docstring
    const reader::T
    const record::SAMRecord
    const offset_first_record::Int
    const header::SAMHeader
    # Number of lines in the header, used to compute line numbers of records
    const n_header_lines::Int
    const copy::Bool
    # Number of records read since the first record, including a failed one
    n_records::Int
    error::Union{Nothing, XAMError}
end

"""
    XAM2.try_init(::Type{SAMReader}, io::Union{IO, AbstractBufReader}; copy::Bool = true)::Union{SAMReader, XAMError}
    XAM2.try_init(::Type{BAMReader}, io::Union{IO, AbstractBufReader}; copy::Bool = true, n_workers::Int = Threads.nthreads())::Union{BAMReader, XAMError}
    XAM2.try_init(::Type{BAMReader}, bgzf::BGZFReader; copy::Bool = true)::Union{BAMReader, XAMError}

Construct a reader like its constructor does, but return an `XAMError` instead of
throwing if the header is malformed.

This function does not throw if the header is malformed. Errors from the underlying IO
are thrown, including errors from `BGZFLib` if the BGZF compression of a BAM file
is malformed. If construction fails, either by returning an error or throwing,
`io` is closed, since the reader takes ownership of it.

# Examples
```julia
reader = XAM2.try_init(SAMReader, open("file.sam"))
if reader isa XAMError
    println("Malformed header: ", sprint(showerror, reader))
else
    # use reader
end
```
"""
function try_init end

function try_init(::Type{SAMReader}, reader::AbstractBufReader; copy::Bool = true)::Union{SAMReader, XAMError}
    header = close_on_throw(() -> read_initial_sam_header!(reader), reader)
    if header isa XAMError
        close(reader)
        return header
    end
    data = MemoryView(header)
    n_header_lines = count(==(UInt8('\n')), data)
    return SAMReader{typeof(reader)}(reader, SAMRecord(), length(data), header, n_header_lines, copy, 0, nothing)
end

function try_init(::Type{SAMReader}, io::IO; copy::Bool = true)::Union{SAMReader, XAMError}
    return try_init(SAMReader, BufReader(io); copy)
end

function SAMReader(io::Union{IO, AbstractBufReader}; copy::Bool = true)
    reader = try_init(SAMReader, io; copy)
    reader isa XAMError && throw(reader)
    return reader
end

function SAMReader(f::Function, args...; kwargs...)
    reader = SAMReader(args...; kwargs...)
    return try
        f(reader)
    finally
        close(reader)
    end
end

function Base.close(reader::SAMReader)
    reader.error = XAMError(XAMErrors.closed_stream)
    return close(reader.reader)
end

get_header(reader::SAMReader) = reader.header

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

function read_initial_sam_header!(reader::AbstractBufReader)::Union{SAMHeader, XAMError}
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
    return try_parse_header(ImmutableMemoryView(get_unflushed(buffer)))
end

function seek_first_record(reader::SAMReader)
    seek(reader.reader, reader.offset_first_record)
    reset!(reader.record)
    reader.n_records = 0
    clear_error!(reader)
    return reader
end

Base.IteratorSize(::Type{<:SAMReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:SAMReader}) = SAMRecord

function Base.iterate(reader::SAMReader, _::Nothing = nothing)
    result = try_read!(reader)
    result === nothing && return nothing
    isa(result, XAMError) && throw(result)
    record = reader.copy ? copy(result) : result
    return (record, nothing)
end

"""
    XAM2.try_read!(reader::SAMReader)::Union{SAMRecord, Nothing, XAMError}

Attempt to read the next record of `reader`, similar to iterating the reader.
If the reader is EOF, return `nothing`.
If the next record is malformed, return an `XAMError`, with the record number,
line and column of the error.
Once this function returns an `XAMError`, subsequent calls to this function
will return the same error until the reader is closed or seeked.

The record is parsed with [`XAM2.try_parse!`](@ref), see its documentation for
what is validated.

This function does not throw if the reader's data is malformed.
"""
function try_read!(reader::SAMReader)::Union{SAMRecord, Nothing, XAMError}
    isnothing(reader.error) || return reader.error
    # Get a view of the current line, stripped from newline
    line_view = get_next_line_view(reader.reader)
    line_view === nothing && return nothing
    reader.n_records += 1
    result = try_parse!(reader.record, _chomp(line_view))
    @inbounds consume(reader.reader, length(line_view))
    # On failure, `try_parse!` has already reset the record
    result isa XAMError && return fail_read!(reader, result)
    return result
end

@noinline function fail_read!(reader::SAMReader, err::XAMError)::XAMError
    record = reader.n_records
    err = XAMError(err.kind; record, line = reader.n_header_lines + record, err.column)
    reader.error = err
    return err
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
