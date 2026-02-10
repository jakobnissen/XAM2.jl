"""
    mutable struct BAMRecord{R}

A single alignment record in a BAM file.

The type parameter `R` is a concrete subtype of `BAMReader` or `BAMWriter`,
and represents the owner of the record, which is needed to retrieve information
on references stored in the record.

`BAMRecord`s are semi-lazily parsed, and do not e.g. validate their CIGAR
before access.

    BAMRecord(io::Union{BAMReader, BAMWriter})::BAMRecord
    
Create an unmapped, empty (but valid) record referring `io`.
"""
mutable struct BAMRecord{R}
    # R can be reader or writer.

    # This data contains, in order, with no unused bytes at the end:
    # 1. Read name including trailing null byte
    # 2. Cigar ops encoded in 32 bits each
    # 3. Seq encoded in 4 bits each
    # 4. Quality
    # 5. Auxiliary data
    const data::Vector{UInt8}
    # We need a ref here because the refids point back to the header,
    # stored in the reader/writer
    const reader_or_writer::R

    # This field is cached. Uninit means: Not computed. Nothing means no cigar.
    # Otherwise, a cigar or invalid data.
    cigar::Union{Uninit, Nothing, BAMCIGAR, CIGARError}

    # fixed-length fields (see BAM specs for the details)
    # Length of alignment record, excluding this 32-bit field
    block_size::UInt32
    # Zero-based index of reference in , -1 if missing
    refID::Int32
    # Zero-based position of leftmost base of alignment, or -1 if missing
    pos::Int32
    # Length of read name plus one (this value includes the trailing null byte)
    l_read_name::UInt8
    # 0xff indicates missing
    mapq::UInt8
    # Bin index (used when indexing BAM files)
    bin::UInt16
    # Number of CIGAR ops - but can be missing, see `is_missing_cigar_field`
    n_cigar_op::UInt16
    flags::SAMFlags
    l_seq::UInt32
    # Same as refid
    next_refID::Int32
    next_pos::Int32
    # 0 if missing, else TLEN as defined in the SAM spec
    tlen::Int32

    global function new_uninitialized_record(
            reader_or_writer::R,
            data::Vector{UInt8}
        ) where {R}
        return new{R}(data, reader_or_writer, Uninit())
    end
end

const BAM_RECORD_CONTENT = [
    0x22, 0x00, 0x00, 0x00, # size of an empty BAM record
    0xff, 0xff, 0xff, 0xff, # unmapped refid
    0xff, 0xff, 0xff, 0xff, # unmapped pos
    0x02, # read name len
    0xff, # mapq unavailable
    0x48, 0x12, # bin 4680
    0x00, 0x00, # no CIGAR ops
    0x04, 0x00, # flag: unmapped,
    0x00, 0x00, 0x00, 0x00, # seq len
    0xff, 0xff, 0xff, 0xff, # unmapped next refid
    0xff, 0xff, 0xff, 0xff, # unmapped next pos
    0x00, 0x00, 0x00, 0x00, # unmapped tlen
] |> to_memory

const BAM_FIXED_BYTES = length(BAM_RECORD_CONTENT)

# Short form: Mapped, read name, strand
function Base.show(io::IO, record::BAMRecord)
    mapped = is_mapped(record) ? "Mapped" : "Unmapped"
    segment = is_first(record) ? "First" : (is_last(record) ? "Last" : "Md")
    qname = query_name(record)
    qn = if qname === nothing
        ""
    else
        ", " * repr(qname)
    end
    return print(io, "BAMRecord<", mapped, ", ", segment, qn, '>')
end

function Base.show(io::IO, ::MIME"text/plain", record::BAMRecord)
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
    ref = reference(record)
    if ref === nothing
        println(buf, "< missing >")
    else
        println(buf, rtruncate(string(ref), width - 12))
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

"""
    mutable struct BAMReader{T <: AbstractBufReader}
    
Lazy, multithreaded reader of BAM files.

`rdr = BAMReader` is an iterator over `BAMRecord{typeof(rdr)}`.
The header can be obtained with [`get_header`](@ref)

    BAMReader(
        io::Union{IO, AbstractBufReader};
        copy::Bool = true,
        n_workers::Int = Threads.nthreads()
    )::BAMReader

    BAMReader(bgzf::BGZFReader; copy::Bool=true)

Construct a `BAMReader` from `io`, by wrapping `io` in a `BGZFLib.BGZFReader`.
This reader takes control of `io`. Mutating `io` after the reader is constructed
is not permitted. Closing the resulting reader will close `io`.

If `copy` is `false`, every record returned by iterating the reader or from
`XAM2.try_read!` will be the same record, overwritten in-place. This improves
efficiency, but new records will invalidate every previous one.
To prevent aliasing, `copy` can be called manually on records.

The number of decompressing worker tasks is set by `n_workers`, which must be at least 1.
"""
mutable struct BAMReader{T <: AbstractBufReader}
    const reader::BGZFReader{T}
    # Lens and names in their own fields, because these are accessed far more often
    # than any other metadata
    const ref_lens::Vector{Int32}
    const ref_names::StringVec
    const header::SAMHeader
    const offset_first_record::VirtualOffset
    const copy::Bool
    error::Union{Nothing, XAMError}

    # These must be the last field so we can use a hack to initialize
    # the record and the reader at the same time.
    # That's also why it can't be const
    record::BAMRecord{BAMReader{T}}

    global function new_bam_reader(reader::BGZFReader{T}, copy::Bool) where {T}
        buffer = UInt8[]
        # TODO: Handle potential errors - is maybe not SAMHeader, can be error
        # also, can throw
        (; header, stringvec, reflens) = parse_bam_header!(reader, buffer)
        vo_at_first_record = virtual_position(reader)

        result = new{T}(
            reader,
            reflens,
            stringvec,
            header,
            vo_at_first_record,
            copy,
            nothing
        )
        result.record = new_bam_record(result, buffer)
        return result
    end
end

# TODO: For both reader and writer
# TODO: Close
# TODO: open(do) method

get_header(reader::BAMReader) = reader.header

const FIXED_BAM_OFFSET = fieldoffset(BAMRecord{BAMReader{CursorReader}}, 4)

"""
    reset!(::BAMRecord)

Resets the `BAMRecord` to its default state, deleting all information from it.
After this function, the record is equal to an empty record created with
`BAMRecord(reader_or_writer)`.
"""
function reset!(record::BAMRecord)
    resize!(record.data, 2)
    @inbounds record.data[1] = UInt8('*')
    @inbounds record.data[2] = 0x00
    record.cigar = Uninit()
    return reset_fixed_fields!(record)
end

function Base.iterate(reader::BAMReader, _::Nothing = nothing)
    result = try_read!(reader)
    result === nothing && return nothing
    isa(result, XAMError) && throw(result)
    record = reader.copy ? copy(result) : result
    return (record, nothing)
end

@inline function reset_fixed_fields!(record::BAMRecord)
    GC.@preserve BAM_RECORD_CONTENT record begin
        dst = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            pointer(BAM_RECORD_CONTENT)::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return record
end

function Base.copy(record::BAMRecord)
    vector_copy = copy(record.data)
    cp = new_uninitialized_record(record.reader_or_writer, vector_copy)
    GC.@preserve cp record begin
        dst = Ptr{UInt8}(pointer_from_objref(cp)) + FIXED_BAM_OFFSET
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            src::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    record.cigar = copy_cigar!(record.cigar, record.data, vector_copy)
    return cp
end

function seek_first_record(reader::BAMReader)
    # TODO: Check code base for when we reset reader.error
    virtual_seek(reader.reader, reader.offset_first_record)
    reset!(reader.record)
    return reader
end

function parse_bam_header!(
        reader::BGZFReader,
        buffer::Vector{UInt8}
    )::Union{
        XAMError, @NamedTuple{
            header::SAMHeader,
            stringvec::StringVec,
            reflens::Vector{Int32},
        }
    }

    memview = fill_bytes_to_data!(reader, buffer, 8)
    isnothing(memview) && return XAMErrors.header_truncated

    magic = unsafe_load_integer_le(UInt32, memview, 1)
    magic == 0x014d4142 || return XAMErrors.header_bad_magic_bytes

    l_text = unsafe_load_integer_le(UInt32, memview, 5)
    l_text > typemax(Int32) % UInt32 && return XAMErrors.header_top_bit_set

    memview = fill_bytes_to_data!(reader, buffer, l_text % Int + 4)
    isnothing(memview) && return XAMErrors.header_truncated

    # Copy the header because we need to store the data in the reader
    header_view = copy(@inbounds(memview[1:(end - 4)]))

    # TODO: This can throw. Is that alright? How do we handle it?
    header = parse_header(header_view)
    stringvec = StringVec()
    reflens = Int32[]

    n_ref = unsafe_load_integer_le(UInt32, memview, l_text % Int + 1)
    n_ref > typemax(Int32) % UInt32 && return XAMErrors.header_top_bit_set

    for _ in 1:(n_ref % Int)
        memview = fill_bytes_to_data!(reader, buffer, 4)
        isnothing(memview) && return XAMErrors.header_truncated

        l_name = unsafe_load_integer_le(UInt32, memview, 1)

        memview = fill_bytes_to_data!(reader, buffer, l_name + 4)
        isnothing(memview) && return XAMErrors.header_truncated

        # Skip null byte
        push_codeunits!(stringvec, memview[1:(l_name - 1)])

        l_ref = unsafe_load_integer_le(Int32, memview, l_name + 1)
        l_ref < 1 && return XAMErrors.header_empty_ref

        push!(reflens, l_ref)
    end

    return (; header, stringvec, reflens)
end

function fill_bytes_to_data!(reader::BAMReader, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, Nothing}
    return fill_bytes_to_data!(reader.reader, reader.record.data, nbytes)
end

function fill_bytes_to_data!(reader::AbstractBufReader, data::Vector{UInt8}, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, Nothing}
    length(data) < nbytes && resize!(data, nbytes)
    remaining = MemoryView(data)[1:nbytes]
    while !isempty(remaining)
        buffer = @something get_nonempty_buffer(reader) return nothing
        mn = min(length(buffer), length(remaining))
        @assert mn > 0
        @assert mn <= length(buffer)
        slice = @inbounds buffer[1:mn]
        copyto!(remaining, slice)
        consume(reader, mn)
        remaining = @inbounds remaining[(mn + 1):end]
    end
    return @inbounds ImmutableMemoryView(data)[1:nbytes]
end

Base.isopen(x::BAMReader) = isopen(x.reader)

function BAMReader(reader::BGZFReader{T}; copy::Bool = true) where {T}
    return new_bam_reader(reader, copy)
end

function BAMReader(io::Union{IO, AbstractBufReader}; copy::Bool = true, n_workers::Int = Threads.nthreads())
    reader = BGZFReader(io; n_workers)
    return new_bam_reader(reader, copy)
end

function Base.show(io::IO, reader::BAMReader)
    print(io, '<')
    summary(io, reader)
    return print(io, '>')
end

function Base.show(io::IO, ::MIME"text/plain", reader::BAMReader)
    summary(io, reader)
    print(io, ':')
    return print(io, "TODO!")
end

function Base.close(reader::BAMReader)
    reader.error = XAMErrors.closed_stream
    return close(reader.reader)
end

BAMRecord(reader::BAMReader) = new_bam_record(reader, UInt8[])

function new_bam_record(reader::BAMReader{T}, v::Vector{UInt8}) where {T}
    resize!(v, 2)
    @inbounds v[1] = UInt8('*')
    @inbounds v[2] = 0x00
    record = new_uninitialized_record(reader, v)
    return reset_fixed_fields!(record)
end

function query_name(record::BAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    view = view[1:Int(record.l_read_name - 0x01)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

function reference_name(record::BAMRecord)::Union{Nothing, AbstractString}
    record.refID < 0 && return nothing
    id = record.refID + 1
    return record.reader_or_writer.ref_names[id]
end

function next_reference_name(record::BAMRecord)::Union{Nothing, AbstractString}
    record.next_refID < 0 && return nothing
    id = record.next_refID + 1
    return record.reader_or_writer.ref_names[id]
end

function reference(record::BAMRecord)
    record.refID < 0 && return nothing
    id = record.refID + 1
    reader = record.reader_or_writer
    # TODO: A single boundscheck here?
    return Reference(
        reader.ref_names[id],
        reader.ref_lens[id]
    )
end

function reference_length(record::BAMRecord)::Union{Int, Nothing}
    record.refID < 0 && return nothing
    id = record.refID + 1
    reader = record.reader_or_writer
    return reader.ref_lens[id] % Int
end

function first_eq_seq_pos(record::BAMRecord)
    iszero(record.l_seq) && return nothing
    n_seq_whole_bytes = record.l_seq >> 1
    offset = record.l_read_name + 4 * record.n_cigar_op
    has_half_byte = isodd(record.l_seq)
    span = (offset + 1):(offset + n_seq_whole_bytes + has_half_byte)
    mem = ImmutableMemoryView(record.data)[span]
    if has_half_byte
        (last_byte, mem) = @inbounds split_last(mem)
        iszero(last_byte & 0xf0) && return record.l_seq % Int
    end
    i = 1
    chunk_size = 128
    @inbounds while i <= length(mem) - chunk_size + 1
        any = false
        for j in i:(i + chunk_size - 1)
            any |= iszero(mem[j] & 0x0f) | iszero(mem[j] & 0xf0)
        end
        any && break
        i += chunk_size
    end
    @inbounds while i <= length(mem)
        (iszero(mem[i] & 0x0f) | iszero(mem[i] & 0xf0)) && return i
        i += 1
    end
    return nothing
end

"""
    sequence(rec::BAMRecord)::Union{Nothing, AbstractString}

Get the query sequence of `rec` as an `AbstractString`, or
`nothing` if the sequence is not available.

The returned `AbstractString`'s type is an implementation detail, but it is guaranteed to:
* Be constructed cheaply in constant time
* Backed by a view into the existing data of the record (so not separately allocated)
* Contain at least 1 codeunit, and having the codeunit type `UInt8`.
* Containing only the `Char`s `"=ACMGRSVTWYHKDBN"` encoded in ASCII
"""
function sequence(record::BAMRecord)::Union{BAMSeq, Nothing}
    iszero(record.l_seq) && return nothing
    offset = record.l_read_name + record.n_cigar_op * 4
    seq_bytes = (record.l_seq + 1) >> 1
    span = (offset + 1):(offset + seq_bytes)
    return BAMSeq(ImmutableMemoryView(record.data)[span], record.l_seq % Int)
end

"""
    quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}

Return the quality as a PHRED-encoded memory view, or `nothing` if the quality is
missing (set to all `0xFF`)

Note that the returned quality does not have a PHRED offset, so e.g. a PHRED
score of 31 is encoded as `UInt8(31)` as opposed to e.g. `UInt8(31 + 33)`.
"""
function quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    seq_bytes = (record.l_seq + 1) >> 1
    offset = record.l_read_name + 4 * record.n_cigar_op + seq_bytes
    qual = ImmutableMemoryView(record.data)[(offset + 1):(offset + record.l_seq)]

    # Missing quality is encoded as all-0xFF
    return if !isempty(qual) && @inbounds first(qual) == 0xff
        nothing
    else
        qual
    end
end

"""
    aux_data(record::BAMRecord)::XAMAuxData.BAM.MutableAuxiliary

Return the auxiliary data as a mutable `XAMAuxData.BAM.Auxiliary`.
If the record does not contain any auxiliary fields, the returned object
is empty.

See the documentation of `XAMAuxData.BAM.Auxiliary` for more details.
"""
function aux_data(record::BAMRecord)
    seq = record.l_seq
    n_seq_bytes = (seq + 1) >> 1
    offset = record.l_read_name + seq + n_seq_bytes + 4 * record.n_cigar_op
    return XAMAuxData.BAM.Auxiliary(record.data, offset + 1)
end

const CIGAR_LUT = let
    u = UInt64(0)
    for (i, byte) in enumerate(b"MIDNSHP=X")
        u |= (byte % UInt64) << (7 * i - 7)
    end
    u
end

function left_position(record::BAMRecord)::Union{Int, Nothing}
    result = record.pos + 1
    return iszero(result) ? nothing : result
end

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError}
    record.cigar === Uninit() || return record.cigar
    encoding = get_cigar_encoding(record)
    if encoding === nothing
        record.cigar = nothing
        return nothing
    end
    return record.cigar = CIGARStrings.try_parse(BAMCIGAR, encoding)
end

@inline function is_missing_encoding(record::BAMRecord, cigar_encoding::ImmutableMemoryView{UInt8})
    length(cigar_encoding) == 8 || return false
    u1 = unsafe_load_integer_le(UInt32, cigar_encoding, 1)
    u2 = unsafe_load_integer_le(UInt32, cigar_encoding, 5)
    return if ((u1 & 0x0f) != 0x04) | ((u2 & 0x0f) != 0x03)
        return false
    else
        ref_length = record.reader_or_writer.ref_lens[record.refID + 1]
        ((u1 >> 4) == record.l_seq) & ((u2 >> 4) == ref_length)
    end
end

function get_cigar_encoding(x::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    iszero(x.n_cigar_op) && return nothing
    lastindex = x.l_read_name + 4 * x.n_cigar_op
    @boundscheck checkbounds(x.data, lastindex)
    view = ImmutableMemoryView(x.data)
    vw = view[(x.l_read_name + 1):lastindex]
    return is_missing_encoding(x, vw) ? load_aux_cigar(x) : vw
end

# Outlined, because it's expected to be rare
@noinline function load_aux_cigar(x::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    v = get(aux_data(x), "CG", nothing)
    isnothing(v) && return v
    return if isa(v, AbstractVector{UInt32})
        parent(v)
    else
        nothing
    end
end

function unsafe_load_fixed_fields(record::BAMRecord, buffer::ImmutableMemoryView{UInt8})
    GC.@preserve record buffer begin
        dst = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            pointer(buffer)::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return nothing
end

"""
    XAM2.try_read(reader::BAMReader)::Union{BAMRecord, Nothing, XAMError}

Attempt to read the next record of `reader`, similar to iterating the reader.
If the reader is EOF, return `nothing`.
If the next record is malformed, return a `XAMError`.
Once this function returns a `XAMError`, subsequent calls to this function
will return the same error until the reader is closed or seeked.

This function does not throw if the reader's data is malformed.
"""
function try_read!(x::BAMReader)::Union{BAMRecord, Nothing, XAMError}
    isnothing(x.error) || return x.error

    # Check if EOF and return that error if so
    buffer = get_nonempty_buffer(x.reader)
    isnothing(buffer) && return nothing

    record = x.record
    vector = record.data

    record.cigar = Uninit()

    # Make room for the fixed bytes
    resize!(vector, BAM_FIXED_BYTES)
    remaining = MemoryView(vector)

    memview = fill_bytes_to_data!(x, BAM_FIXED_BYTES)
    if isnothing(memview)
        x.error = XAMErrors.truncated_file
        return x.error
    end

    unsafe_load_fixed_fields(record, memview)

    # Read all the variable length.
    # The block size is the total length, except the first 4 byte field,
    # account for that.
    remaining_bytes = record.block_size - BAM_FIXED_BYTES + 4

    # Ensure the total field length is not shorter than the sum of the
    # variable field names
    seq_bytes = (record.l_seq + UInt32(1)) >> 1
    minimum_bytes = record.l_read_name + record.l_seq + seq_bytes + 4 * record.n_cigar_op

    if remaining_bytes < minimum_bytes
        x.error = XAMErrors.invalid_block_length
        return x.error
    end

    # Read all the remaining bytes, overwriting the fixed field bytes,
    # which we have already moved into the integer fields of the struct
    resize!(vector, remaining_bytes)
    remaining = MemoryView(vector)

    memview = fill_bytes_to_data!(x, remaining_bytes)
    if isnothing(memview)
        x.error = XAMErrors.truncated_file
        return x.error
    end
    return x.record
end

# TODO: Figure out how this should be constructed
mutable struct BAMWriter{T <: AbstractBufWriter}
    const writer::BGZFWriter{T}
    const ref_lens::Vector{Int32}
    const ref_names::StringVec
    const offset_of_refname::Dict{StringView{ImmutableMemoryView{UInt8}}, UInt32}
    record::BAMRecord{BAMWriter{T}}

    global function new_bam_writer(
            writer::BGZFWriter{T},
            ref_lens::Vector{Int32},
            ref_names::StringVec,
            offset_of_refname::Dict{StringView{ImmutableMemoryView{UInt8}}, UInt32},
        ) where {T}
        w = new{T}(writer, ref_lens, ref_names, offset_of_refname)
        record = new_uninitialized_record(w, UInt8[])
        w.record = reset!(record)
        return w
    end
end

function BAMWriter(
        writer::Union{AbstractBufWriter, IO},
        header::ImmutableMemoryView{UInt8},
        refnames::AbstractVector{<:AbstractString},
        reflens::AbstractVector{<:Integer};
        n_workers::Int = min(4, Threads.nthreads()),
        compress_level::Int = 6,
    )
    bgzf = BGZFWriter(writer; n_workers, compress_level)
    return BAMWriter(bgzf, header, refnames, reflens)
end

# TODO: Reconsider how this should be constructed. Currently, the header
# can conflict with refnames and reflens. That doesn't make sense, is bad API.
function BAMWriter(
        writer::BGZFWriter{T},
        header::ImmutableMemoryView{UInt8},
        refnames::AbstractVector{<:AbstractString},
        reflens::AbstractVector{<:Integer},
    ) where {T}
    if length(header) > typemax(Int32)
        throw(ArgumentError("Header must be 2147483647 bytes or shorter"))
    end
    if length(refnames) != length(reflens)
        throw(ArgumentError("Length of refnames and reflens must be the same"))
    end
    if length(refnames) > typemax(Int32)
        throw(ArgumentError("BAM format supports at most 2147483647 reference sequences"))
    end
    ref_lens = Vector{Int32}(undef, length(reflens))
    for i in eachindex(ref_lens, reflens)
        reflen = reflens[i]
        if reflen < 1 || reflen > typemax(Int32)
            throw(ArgumentError("Reflen must be in 1:2147483647"))
        end
        ref_lens[i] = reflen % Int32
    end

    # Build the dict of refnames
    stringvec = StringVec()
    for refname in refnames
        cu = codeunits(refname)
        # Max refname length includes trailing null byte
        if length(cu) > (typemax(Int32) - 1)
            throw(ArgumentError("Reference name cannot exceed 2147483646 bytes"))
        end
        push_codeunits!(stringvec, cu)
    end
    offset_of_refname = Dict{StringView{ImmutableMemoryView{UInt8}}, UInt32}()
    for i in eachindex(stringvec)
        offset_of_refname[stringvec[i]] = i - 1
    end

    @assert length(stringvec) == length(ref_lens)

    # Write the header:
    # Write magic bytes BAM\1
    write(writer, 0x014d4142)
    write(writer, length(header) % UInt32)
    write(writer, header)
    write(writer, length(refnames) % UInt32)

    # Write reference sequences
    for (refname, ref_len) in zip(stringvec, ref_lens)
        cu = codeunits(refname)
        write(writer, (length(cu) + 1) % UInt32)
        write(writer, cu, 0x00)
        write(writer, ref_len % UInt32)
    end

    return new_bam_writer(writer, ref_lens, stringvec, offset_of_refname)
end

function Base.show(io::IO, writer::BAMWriter)
    print(io, '<')
    summary(io, writer)
    return print(io, '>')
end

function Base.show(io::IO, ::MIME"text/plain", writer::BAMWriter)
    summary(io, writer)
    print(io, ':')
    return print(io, "TODO!")
end

function Base.close(writer::BAMWriter)
    return close(writer)
end

function Base.write(writer::BAMWriter, record::BAMRecord)
    # Write fixed bytes
    buffer = get_buffer(writer)
    if length(buffer) < BAM_FIXED_BYTES
        grow_buffer(writer)
        buffer = get_buffer(writer)
        # This is unreachable, because calling `grow_buffer` on a BGZFWriter
        # will always make sure the next buffer is a full BGZF block,
        # which is 64 KiB.
        length(buffer) < BAM_FIXED_BYTES && error("unreachable")
    end
    GC.@preserve record buffer begin
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            pointer(buffer)::Ptr{Cchar},
            src::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    @inbounds consume(writer.writer, BAM_FIXED_BYTES)
    # Write variable bytes
    return BAM_FIXED_BYTES + write(writer.writer, record.data)
end

function store_dummy_bamcigar!(mem::MutableMemoryView{UInt8}, l_seq::UInt32, l_ref::UInt32)
    max(l_seq, l_ref) > 0x0fff_ffff && error() # TODO
    length(mem) == 8 || error() # TODO
    GC.@preserve mem begin
        ptr = Ptr{UInt32}(pointer(mem))
        unsafe_store!(ptr, htol((l_seq << 4 | 0x04)))
        unsafe_store!(ptr + 4, htol((l_ref << 4 | 0x03)))
    end
    return nothing
end

function get_index_of_refname(record::BAMRecord{<:BAMWriter}, name::Union{AbstractString, Nothing})
    name === nothing && return Int32(-1)
    id = get(record.reader_or_writer.offset_of_refname, name, nothing)
    return if id === nothing
        error() # TODO: Throw a XAM error
    else
        Int32(id)
    end
end

function Base.copy!(dst::BAMRecord{<:BAMWriter}, src::SAMRecord)
    # Copy over basic fields
    dst.pos = (src.pos - 1) % Int32
    dst.next_pos = (src.pnext - 1) % Int32
    dst.mapq = src.mapq
    dst.flags = src.flags
    dst.tlen = src.tlen

    # Reference IDs
    dst.refID = get_index_of_refname(dst, reference_name(src))
    dst.next_refID = get_index_of_refname(dst, next_reference_name(src))

    # Mapping span and bin
    dst.bin = let
        ms = mapping_span(src)
        if ms === nothing
            UInt16(4680)
        else
            span = (first(ms) % Int32):(last(ms) % Int32)
            reg2bin(span)
        end
    end

    # Set dst.n_cigar_op, handling if there is more than typemax(UInt16) ops,
    # in which case specs say to store a particularly dummy value of 2 ops,
    # then store it in the aux. Here, we simply store if the cigar is too big,
    # and dst.n_cigar_op
    cig = cigar(src)
    cigar_too_big = false
    dst.n_cigar_op = if cig === nothing
        0
    elseif cig isa XAMError
        throw(XAMErrors.invalid_cigar) # TODO: Real error <: Exception type
    else
        cigar_too_big = length(cig) > typemax(UInt16)
        # If CIGAR is to big, we store this dummy value and store the real cigar in an AUX tag
        cigar_too_big ? UInt16(2) : length(cig) % UInt16
    end

    # Sequence length
    seq = sequence(src)
    dst.l_seq = seq === nothing ? 0 : length(seq) % UInt32

    # Resize dst's vector to hold all variable-length fields, except auxdata, which we add later
    v_len = src.qname_len + 1 + 4 * dst.n_cigar_op + ((dst.l_seq + 1) >>> 1) + dst.l_seq
    resize!(dst.data, v_len)
    mem = MemoryView(dst.data)

    # Copy over qname + trailing zero byte
    qname = query_name(src)
    if qname === nothing
        mem[1] = UInt8('*')
        mem[2] = 0x00
        dst.l_read_name = 0x02
    else
        cu = codeunits(qname)
        copyto!(mem, cu)
        mem[length(cu) + 1] = 0x00
        dst.l_read_name = UInt8(length(cu) + 1)
    end
    offset = dst.l_read_name

    # Copy over cigar, either the dummy value if it's too big, or else the actuial cihar
    if cig isa CIGAR
        if cigar_too_big
            # If c isa CIGAR, the record must be mapped, and so refID cannot be negative.
            reflen = dst.reader_or_writer.ref_lens[dst.refID + 1]
            store_dummy_bamcigar!(mem[(offset + 1):(offset + 8)], dst.l_seq, reflen % UInt32)
            dst.cigar = Uninit()
            offset += 8
        else
            dst_mem = mem[(offset + 1):(offset + 4 * dst.n_cigar_op)]
            dst.cigar = CIGARStrings.encode!(dst_mem, BAMCIGAR, cig)
            offset += 4 * dst.n_cigar_op
        end
    end

    # Copy over sequence
    if seq !== nothing
        bamseq!(mem[(offset + 1):end], codeunits(seq))
        offset += (dst.l_seq + 1) >>> 1
    end

    # Copy over quality
    qual = quality_scores(src)
    if qual !== nothing
        seq === nothing && error() # TODO
        v = mem[(offset + 1):end]
        cu = codeunits(qual)
        for i in eachindex(cu, v)
            v[i] = cu[i] - UInt8(33)
        end
    end

    # Copy over aux, if any
    if src.aux_offset < length(src.data)
        dst_aux = XAMAuxData.BAM.Auxiliary(dst.data, length(dst.data) + 1)
        copy!(dst_aux, aux_data(src))
    end

    # If cigar was too big, we copy over its binary encoding to the AUX directly here.
    if cigar_too_big
        cig isa CIGAR || throw(AssertionError())
        aux = aux_data(dst)
        if haskey(aux, AuxTag("CG"))
            error() # TODO
        end

        # Manually write the cigar in, without allocating an array
        cigar_len = 4 * length(cig) + 6
        oldlen = length(dst.data)
        resize!(dst.data, oldlen + cigar_len)
        mem = MemoryView(dst.data)[(oldlen + 1):end]
        mem[1] = UInt8('B')
        mem[2] = UInt8('i')
        GC.@preserve mem begin
            ptr = Ptr{UInt32}(pointer(mem)) + 2
            unsafe_store!(ptr, UInt32(cigar_len))
            for element in c
                u = reinterpret(UInt8, element.op) % UInt32
                u |= (element.len % UInt32) << 4
                ptr += 4
                unsafe_store!(ptr, u)
            end
        end
    end

    # Finally, the total block size can be computed
    dst.block_size = (length(dst.data) + BAM_FIXED_BYTES - 4) % UInt32

    return dst
end

function append_number_field!(
        data::Vector{UInt8},
        n::Unsigned,
        max_digits::Int
    )::Int
    old_len = length(data)
    resize!(data, old_len + max_digits + 1)
    n = inbounds_write_base10!(@inbounds(MemoryView(data)[(old_len + 1):end]), n)
    resize!(data, old_len + n + 1)
    @inbounds data[end] = UInt8('\t')
    return old_len + n + 1
end

min_string(::Int64) = "-9223372036854775808"
min_string(::Int32) = "-2147483648"

function append_number_field!(
        data::Vector{UInt8},
        n::Union{Int32, Int64},
        max_digits::Int
    )::Int
    if n == typemin(n)
        s = min_string(typemin(n))
        old_len = length(data)
        sizehint!(data, old_len + sizeof(s) + 1; shrink = true)
        append!(data, codeunits(s))
        @inbounds data[end] = UInt8('\t')
    else
        if n < zero(n)
            push!(data, UInt8('-'))
        end
        return append_number_field!(data, unsigned(abs(n)), max_digits) + (n < zero(n))
    end
end

function append_optional_string_field!(data::Vector{UInt8}, field::Union{Nothing, AbstractString})::Int
    cu = if field === nothing
        push!(data, UInt8('*'))
        1
    else
        append!(data, codeunits(field))
        ncodeunits(field)
    end
    push!(data, UInt8('\t'))
    return cu + 1
end

function Base.copy!(dst::SAMRecord, src::BAMRecord)
    data = dst.data
    empty!(data)

    # QNAME
    n = append_optional_string_field!(data, query_name(src))
    dst.qname_len = (n % UInt8) - 0x01

    # FLAG
    append_number_field!(data, src.flags.x, 5)
    dst.flags = src.flags

    # RNAME
    rname_start = (length(data) + 1) % UInt32
    rname_len_p1 = append_optional_string_field!(data, reference_name(src))
    dst.rname_span = (rname_start, (rname_start + rname_len_p1 - 1) % UInt32)

    # POS
    dst.pos = src.pos + UInt32(1)
    append_number_field!(data, dst.pos, 10)

    # MAPQ
    append_number_field!(data, src.mapq, 10)
    dst.mapq = src.mapq

    # CIGAR
    dst.cigar_offset = length(data) % UInt32
    cig = cigar(src)
    if cig isa CIGARError
        throw(XAMErrors.invalid_cigar)
    elseif cig === nothing
        push!(data, UInt8('*'))
    else
        dst.cigar = CIGARStrings.encode_append!(data, CIGAR, cig)
    end
    push!(data, UInt8('\t'))

    # RNEXT
    if src.refID == src.next_refID
        # If RNAME and RNEXT is the same, we write it as `=`
        idx = length(data + 1) % UInt32
        dst.rnext_span = (idx, idx)
        push!(data, UInt8("="))
        push!(data, UInt8("\t"))
    else
        rnext_start = (length(data) + 1) % UInt32
        rnext_name_len_p1 = append_optional_string_field!(data, next_reference_name(src))
        dst.rnext_span = (rnext_start, (rnext_start + rnext_name_len_p1 - 1) % UInt32)
    end

    # PNEXT
    dst.pnext = src.next_pos + UInt32(1)
    append_number_field!(data, dst.pnext, 10)

    # TLEN
    dst.tlen = src.tlen
    append_number_field!(data, src.tlen, 10)

    # SEQ
    dst.seq_offset = length(data) % UInt32
    append_optional_string_field!(data, sequence(src))

    # QUAL
    dst.qual_offset = length(data) % UInt32
    qual = quality_scores(src)
    if qual === nothing
        push!(data, UInt8('*'))
    else
        old_len = length(data)
        resize!(data, old_len + length(qual))
        for (i, byte) in enumerate(qual)
            data[old_len + i] = byte + UInt8(33)
        end
    end

    # Auxiliary fields
    dst.aux_offset = (length(data) + 1) % UInt32 # compensate for tab
    src_aux = aux_data(src)
    if !isempty(src_aux)
        push!(data, UInt8('\t'))
        dst_aux = XAMAuxData.SAM.Auxiliary(data, length(data) + 1)
        copy!(dst_aux, src_aux)
    end

    return dst
end
