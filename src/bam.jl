"""
    mutable struct BAMRecord

A single alignment record in a BAM file.

Every record refers to a [`References`](@ref), which is used to look up the
names and lengths of the references that the record's reference IDs point to.
Records read from a `BAMReader` share the reader's `References`.

`BAMRecord`s are semi-lazily parsed: When reading a record, its structure is validated,
e.g. that its reference IDs and positions are in bounds, and that its fields fit in
the record. The content of variable-length fields, like the CIGAR, is validated
when accessed. Accessors never throw on malformed data: Malformed content is returned
as an `XAMError`, e.g. by [`cigar`](@ref) for a malformed CIGAR.

Auxiliary data is returned as an `XAMAuxData.BAM.Auxiliary` by [`aux_data`](@ref),
which returns an error if it is malformed. [`aux_data_unchecked`](@ref) skips this check,
and the returned object may then throw when accessed.

    BAMRecord(refs::References)::BAMRecord
    BAMRecord(io::Union{BAMReader, BAMWriter})::BAMRecord

Create an unmapped, empty (but valid) record referring to `refs`, or to the
references of `io`.
"""
mutable struct BAMRecord
    # Structural invariants. Every `BAMRecord` satisfies these, so accessors may rely on them,
    # e.g. using `@inbounds`. Records from `try_read!` are checked by `check_fixed_fields`,
    # and for the null terminator of the read name. Records created any other way are valid
    # by construction. The fields of `BAMRecord` are internal, so users can't violate these.
    # 1. `length(data) <= typemax(UInt32) - 32`, so the block size fits in a `UInt32`
    #    (see `block_size`)
    # 2. The variable-length fields fit in `data`, i.e.
    #    `l_read_name + 4 * n_cigar_op + cld(l_seq, 2) + l_seq <= length(data)`,
    #    computed without overflow. The rest of `data` is auxiliary data.
    # 3. `l_read_name >= 1`, and `data[l_read_name] == 0x00`
    # 4. `refID` and `next_refID` are in `-1:(length(refs) - 1)`
    # 5. `pos` and `next_pos` are in `-1:(typemax(Int32) - 1)`. The upper bound is the
    #    largest 0-based position that SAM's 1-based POS and PNEXT can represent.
    # The content of the variable-length fields (CIGAR, quality, auxiliary data) is not
    # validated when reading, but when accessed. Accessors return malformed content
    # as an error value, and never throw.

    # This data contains, in order, with no unused bytes at the end:
    # 1. Read name including trailing null byte
    # 2. Cigar ops encoded in 32 bits each
    # 3. Seq encoded in 4 bits each
    # 4. Quality
    # 5. Auxiliary data
    const data::Vector{UInt8}
    # The reference IDs of the record index into these references
    const refs::References

    # Fixed-length fields (see BAM specs for the details), in the same layout as in BAM,
    # such that they can be copied with `memcpy`.
    # The block size field that precedes them in BAM is not stored, see `block_size`.
    # Zero-based index of reference in , -1 if missing
    refID::Int32
    # Zero-based position of leftmost base of alignment, or -1 if missing
    pos::Int32
    # Length of read name plus one (this value includes the trailing null byte)
    l_read_name::UInt8
    # 0xff indicates missing
    mapq::UInt8
    # Bin index, used when indexing BAM files. This is never read: It may be stale, e.g.
    # if the record was edited, so the writer computes the bin instead, see `compute_bin`.
    # The field is kept so the fixed-length fields have the same layout as in BAM.
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

    global function new_uninitialized_record(refs::References, data::Vector{UInt8})
        return new(data, refs)
    end
end

const BAM_RECORD_CONTENT = [
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
] |> Tuple

# Number of bytes of fixed-length fields stored in `BAMRecord`. In BAM, they are
# preceded by the 4-byte block size.
const BAM_FIXED_BYTES = length(BAM_RECORD_CONTENT)

# The BAM block size: The number of bytes of the record after the block size field.
# Invariant: This fits in a `UInt32`.
block_size(record::BAMRecord) = (length(record.data) + BAM_FIXED_BYTES) % UInt32

# Short form: Mapped, read name, strand
function Base.show(io::IO, record::BAMRecord)
    mapped = is_mapped(record) ? "Mapped" : "Unmapped"
    segment = is_first(record) ? "First" : (is_last(record) ? "Last" : "Middle")
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
    elseif cig isa XAMError
        println(buf, "< malformed CIGAR >")
    else
        str = String(MemoryView(CIGAR(cig)))
        println(buf, rtruncate(str, width - 12))
    end

    # Auxiliary data
    aux = aux_data(record)
    if aux isa XAMError
        println(buf, "  Auxiliary: < malformed >")
    elseif !isempty(aux)
        println(buf, "  Auxiliary:")
        auxbuf = IOBuffer()
        show(auxbuf, MIME"text/plain"(), aux)
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

`rdr = BAMReader` is an iterator over `BAMRecord`.
The header can be obtained with [`get_header`](@ref), and the reference
sequences with [`references`](@ref).

    BAMReader(
        io::Union{IO, AbstractBufReader};
        copy::Bool = true,
        n_workers::Int = Threads.nthreads(),
        strict_header::Bool = false,
    )::BAMReader

    BAMReader(bgzf::BGZFReader; copy::Bool = true, strict_header::Bool = false)

Construct a `BAMReader` from `io`, by wrapping `io` in a `BGZFLib.BGZFReader`,
reading and validating the header.
This reader takes control of `io`. Mutating `io` after the reader is constructed
is not permitted. Closing the resulting reader will close `io`.
To obtain the `BGZFReader`, use `BufferIO.take_inner(::BAMReader)`.

If the binary header is malformed, an `XAMError` is thrown.
See [`XAM2.try_init`](@ref) for a version that returns the error instead.
If construction fails, `io` (or `bgzf`) is not closed, and remains owned by the caller.
Its position is then unspecified.
By default, a malformed text header, or one whose `@SQ` lines disagree with the
binary references, does not prevent reading records, see [`get_header`](@ref).
If `strict_header` is `true`, such a header is instead an error when constructing
the reader, so `get_header` is guaranteed to return a `SAMHeader`.
See also [`XAM2.try_validate`](@ref).

Use `seek(reader, offset::VirtualOffset)` to seek to a record, e.g. to the
`offset` of an `XAMError` returned when reading a record. Seeking to an offset outside
the file, or to a block offset outside its block, throws a `BGZFLib.BGZFError`.

Malformed BGZF compression does not throw when constructing the reader or reading
records: it is returned as an `XAMError` of kind `XAMErrors.bgzf_error`. Errors from
the underlying IO are thrown.

If `copy` is `false`, every record returned by iterating the reader or from
`XAM2.try_read!` will be the same record, overwritten in-place. This improves
efficiency, but new records will invalidate every previous one.
To prevent aliasing, `copy` can be called manually on records.

The number of decompressing worker tasks is set by `n_workers`, which must be at least 1.
"""
mutable struct BAMReader{T <: AbstractBufReader}
    const reader::BGZFReader{T}
    const refs::References
    # The text header, or the error from parsing it. It's not needed to read records,
    # since the canonical information needed can be found in the References
    const header::Union{SAMHeader, XAMError}
    # Cached to seeking to start does not need to re-parse the header
    const offset_first_record::VirtualOffset
    const copy::Bool
    const record::BAMRecord
    error::Union{Nothing, XAMError}
end

function try_init(
        ::Type{BAMReader},
        reader::BGZFReader{T};
        copy::Bool = true,
        strict_header::Bool = false,
    )::Union{BAMReader{T}, XAMError} where {T}
    # The reader only takes ownership of `reader` on success, so it's not closed on failure
    buffer = UInt8[]
    parsed = parse_bam_header!(reader, buffer)
    parsed isa XAMError && return parsed
    (; header, refs) = parsed
    strict_header && header isa XAMError && return header
    vo_at_first_record = virtual_position(reader)
    record = new_bam_record(refs, buffer)
    return BAMReader{T}(reader, refs, header, vo_at_first_record, copy, record, nothing)
end

function try_init(
        ::Type{BAMReader},
        io::Union{IO, AbstractBufReader};
        copy::Bool = true,
        n_workers::Int = Threads.nthreads(),
        strict_header::Bool = false,
    )::Union{BAMReader, XAMError}
    # If the constructor throws, it has not spawned any workers, so there is nothing to clean up
    bgzf = BGZFReader(io; n_workers)
    # We own `bgzf`, but the caller owns `io`, so on failure, stop the workers of `bgzf`
    # without closing `io`.
    result = try
        try_init(BAMReader, bgzf; copy, strict_header)
    catch
        BufferIO.take_inner(bgzf)
        rethrow()
    end
    result isa XAMError && BufferIO.take_inner(bgzf)
    return result
end

get_header(reader::BAMReader) = reader.header
references(reader::BAMReader) = reader.refs

const FIXED_BAM_OFFSET = fieldoffset(BAMRecord, 3)

@assert sizeof(BAMRecord) == FIXED_BAM_OFFSET + BAM_FIXED_BYTES

"""
    reset!(record::BAMRecord)::BAMRecord

Reset `record` to its default state, deleting all information from it, and return it.
After this function, the record is equal to an empty record created with
`BAMRecord(record.refs)`. The record keeps its allocated memory, so it can be reused
without allocating.

See also: [`set_record!`](@ref)
"""
function reset!(record::BAMRecord)
    resize!(record.data, 2)
    @inbounds record.data[1] = UInt8('*')
    @inbounds record.data[2] = 0x00
    return reset_fixed_fields!(record)
end

Base.IteratorSize(::Type{<:BAMReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:BAMReader}) = BAMRecord

function Base.iterate(reader::BAMReader, _::Nothing = nothing)
    result = try_read_record!(reader, reader.copy)
    result === nothing && return nothing
    isa(result, XAMError) && throw(result)
    return (result, nothing)
end

@inline function reset_fixed_fields!(record::BAMRecord)
    # Content is stored as a tuple, not as a pointer to a const global `Memory`,
    # because that pointer may be constant folded, which is invalid in a binary
    # compiled with juliac.
    GC.@preserve record begin
        dst = Ptr{typeof(BAM_RECORD_CONTENT)}(pointer_from_objref(record) + FIXED_BAM_OFFSET)
        unsafe_store!(dst, BAM_RECORD_CONTENT)
    end
    return record
end

function Base.copy(record::BAMRecord)
    cp = new_uninitialized_record(record.refs, copy(record.data))
    return copy_fixed_fields!(cp, record)
end

# The fixed-length fields, except `bin`, which is not part of the record's content
function fixed_fields(x::BAMRecord)
    return (
        x.refID, x.pos, x.l_read_name, x.mapq, x.n_cigar_op,
        x.flags, x.l_seq, x.next_refID, x.next_pos, x.tlen,
    )
end

# The references are compared by content, since a record's reference IDs refer to them.
function Base.:(==)(a::BAMRecord, b::BAMRecord)
    fixed_fields(a) == fixed_fields(b) || return false
    a.data == b.data || return false
    return a.refs === b.refs || a.refs == b.refs
end

# The references are not hashed, since hashing them is slow, and it's rarely useful
# to distinguish between records that only differ by their references.
function Base.hash(x::BAMRecord, h::UInt)
    return hash(x.data, hash(fixed_fields(x), h ⊻ (0x7b3d0e5f9a64c812 % UInt)))
end

function copy_fixed_fields!(dst::BAMRecord, src::BAMRecord)
    GC.@preserve dst src begin
        dst_ptr = Ptr{UInt8}(pointer_from_objref(dst)) + FIXED_BAM_OFFSET
        src_ptr = Ptr{UInt8}(pointer_from_objref(src)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            dst_ptr::Ptr{Cchar},
            src_ptr::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return dst
end

function seek_first_record(reader::BAMReader)
    return seek(reader, reader.offset_first_record)
end

function Base.seek(reader::BAMReader, offset::VirtualOffset)
    # After `take_inner`, the `BGZFReader` is still open, but is no longer owned by `reader`
    is_closed(reader) && throw(IOError(IOErrorKinds.ClosedIO))
    # Seeking fills the buffer, which may fail on malformed BGZF data. Then, the error
    # is returned on the next read, as if reading from `offset` failed.
    result = catch_data_error(() -> seek(reader.reader, offset))
    reset!(reader.record)
    clear_error!(reader)
    if result isa BGZFError
        reader.error = XAMError(result; offset)
    end
    return reader
end

# Fill `nbytes` of the header into `buffer`, or return the error.
function fill_header_bytes!(reader::BGZFReader, buffer::Vector{UInt8}, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, XAMError}
    result = fill_bytes_to_data!(reader, buffer, nbytes)
    result === nothing && return XAMError(XAMErrors.header_truncated)
    result isa BGZFError && return XAMError(result)
    return result
end

function parse_bam_header!(
        reader::BGZFReader,
        buffer::Vector{UInt8}
    )::Union{XAMError, @NamedTuple{header::Union{SAMHeader, XAMError}, refs::References}}

    memview = fill_header_bytes!(reader, buffer, 8)
    memview isa XAMError && return memview

    magic = unsafe_load_integer_le(UInt32, memview, 1)
    magic == 0x014d4142 || return XAMError(XAMErrors.header_bad_magic_bytes)

    l_text = unsafe_load_integer_le(UInt32, memview, 5)
    l_text > typemax(Int32) % UInt32 && return XAMError(XAMErrors.header_top_bit_set)

    memview = fill_header_bytes!(reader, buffer, l_text % Int + 4)
    memview isa XAMError && return memview

    # The text may be padded with null bytes, which are not part of the header.
    # Copy the header because we need to store the data in the reader.
    text_len = something(findlast(!iszero, @inbounds(memview[1:(end - 4)])), 0)
    header_view = copy(@inbounds(memview[1:text_len]))

    # Errors in the text header are stored, not returned, since it's not needed to read
    # records.
    header = try_parse_header(header_view)
    stringvec = StringVec()
    reflens = Int32[]

    n_ref = unsafe_load_integer_le(UInt32, memview, l_text % Int + 1)
    n_ref > typemax(Int32) % UInt32 && return XAMError(XAMErrors.header_top_bit_set)

    for _ in 1:(n_ref % Int)
        memview = fill_header_bytes!(reader, buffer, 4)
        memview isa XAMError && return memview

        # Name length includes the trailing null byte
        l_name = unsafe_load_integer_le(UInt32, memview, 1)
        in(l_name, 2:(typemax(Int32) % UInt32)) || return XAMError(XAMErrors.invalid_reference)
        l_name = l_name % Int

        memview = fill_header_bytes!(reader, buffer, l_name + 4)
        memview isa XAMError && return memview

        @inbounds(memview[l_name]) == 0x00 || return XAMError(XAMErrors.invalid_reference)
        # The offsets of `StringVec` are `UInt32`
        if length(stringvec.data) + l_name - 1 > typemax(UInt32)
            return XAMError(XAMErrors.invalid_reference)
        end
        push_codeunits!(stringvec, @inbounds memview[1:(l_name - 1)])

        l_ref = unsafe_load_integer_le(Int32, memview, l_name + 1)
        push!(reflens, l_ref)
    end

    refs = try_build_references(stringvec, reflens)
    refs isa XAMErrorKind && return XAMError(refs)
    # The binary references are used for records, and the text header must not disagree
    if header isa SAMHeader
        header = @something try_validate(header, refs) header
    end
    return @NamedTuple{header::Union{SAMHeader, XAMError}, refs::References}((header, refs))
end

function fill_bytes_to_data!(reader::BAMReader, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, Nothing, BGZFError}
    return fill_bytes_to_data!(reader.reader, reader.record.data, nbytes)
end

# Return the first `nbytes` of `data`, after copying them from `reader`.
# Return `nothing` if `reader` is EOF before that, or a `BGZFError` if the BGZF data
# is malformed.
function fill_bytes_to_data!(reader::AbstractBufReader, data::Vector{UInt8}, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, Nothing, BGZFError}
    length(data) < nbytes && resize!(data, nbytes)
    remaining = MemoryView(data)[1:nbytes]
    while !isempty(remaining)
        buffer = try_get_nonempty_buffer(reader)
        buffer isa BGZFError && return buffer
        buffer === nothing && return nothing
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

# Like `BufferIO.get_nonempty_buffer`, but return a `BGZFError` on malformed BGZF data.
# The buffer is only filled when empty, which is about once per BGZF block, so the cost
# of `try` in `try_fill_buffer` is negligible.
function try_get_nonempty_buffer(reader::AbstractBufReader)::Union{Nothing, ImmutableMemoryView{UInt8}, BGZFError}
    buffer = get_buffer(reader)::ImmutableMemoryView{UInt8}
    isempty(buffer) || return buffer
    result = try_fill_buffer(reader)
    result isa BGZFError && return result
    buffer = get_buffer(reader)::ImmutableMemoryView{UInt8}
    return isempty(buffer) ? nothing : buffer
end

# Fill the buffer, returning malformed BGZF data as a value. See `catch_data_error`.
try_fill_buffer(reader::AbstractBufReader) = catch_data_error(() -> fill_buffer(reader))

# Call `f()`, returning a `BGZFError` from malformed BGZF data as a value, and
# propagating other exceptions, e.g. from the underlying IO.
# After a `BGZFError`, the `BGZFReader` is in an error state until it is seeked.
# TODO: Remove this when BufferIO supports errors as values.
@noinline function catch_data_error(f)
    return try
        f()
    catch e
        (e isa BGZFError && is_data_error(e)) || rethrow()
        e
    end
end

Base.isopen(x::BAMReader) = !is_closed(x) && isopen(x.reader)

function BAMReader(reader::BGZFReader; copy::Bool = true, strict_header::Bool = false)
    result = try_init(BAMReader, reader; copy, strict_header)
    result isa XAMError && throw(result)
    return result
end

function BAMReader(
        io::Union{IO, AbstractBufReader};
        copy::Bool = true,
        n_workers::Int = Threads.nthreads(),
        strict_header::Bool = false,
    )
    result = try_init(BAMReader, io; copy, n_workers, strict_header)
    result isa XAMError && throw(result)
    return result
end


function BAMReader(f::Function, args::Vararg{Any, N}; kwargs...) where {N}
    reader = BAMReader(args...; kwargs...)
    return try
        f(reader)
    finally
        close(reader)
    end
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

# Closing twice is a no-op, so that closing after `take_inner` does not close the
# taken reader, e.g. in the `do`-block form of `BAMReader`.
function Base.close(reader::BAMReader)
    is_closed(reader) && return nothing
    reader.error = XAMError(XAMErrors.closed_stream)
    return close(reader.reader)
end

"""
    virtual_position(reader::BAMReader)::VirtualOffset

Get the virtual offset of the next record to be read from `reader`, or of the end
of the file if `reader` is EOF. Seeking `reader` to this offset with `seek` resumes
reading from that record.
After `reader` has returned an error, the position is unspecified.

If the next record starts at the beginning of a BGZF block, the offset may instead be the
end of the previous block, which is equivalent when seeking, and vice versa.

Throw an `IOError` of kind `ClosedIO` if `reader` is closed.
"""
function BGZFLib.virtual_position(reader::BAMReader)
    # After `take_inner`, the `BGZFReader` is still open, but is no longer owned by `reader`
    is_closed(reader) && throw(IOError(IOErrorKinds.ClosedIO))
    return virtual_position(reader.reader)
end

"""
    BufferIO.take_inner(reader::BAMReader{T})::BGZFReader{T}

Close `reader` and take its `BGZFLib.BGZFReader`, without closing it.
The `BGZFReader` can itself be closed, or passed to
`BufferIO.take_inner` to get the underlying compressed stream.

The returned BGZF reader contains the decompressed BAM data,
and is positioned after the last
record read by `reader`, or at the first record if no records were read.
After `reader` returned an error, its position is unspecified.

Calling this on an already closed `reader` throws a `BufferIO.IOError` of kind `ClosedIO`.
"""
function BufferIO.take_inner(reader::BAMReader)
    is_closed(reader) && throw(IOError(IOErrorKinds.ClosedIO))
    reader.error = XAMError(XAMErrors.closed_stream)
    return reader.reader
end

BAMRecord(refs::References) = new_bam_record(refs, UInt8[])

function new_bam_record(refs::References, v::Vector{UInt8})
    resize!(v, 2)
    @inbounds v[1] = UInt8('*')
    @inbounds v[2] = 0x00
    record = new_uninitialized_record(refs, v)
    return reset_fixed_fields!(record)
end

function query_name(record::BAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    # Invariant: l_read_name >= 1, and the read name fits in the data
    view = @inbounds view[1:(record.l_read_name - 1)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

# Invariant: refID and next_refID are in -1:(length(refs) - 1)
function reference_name(record::BAMRecord)::Union{Nothing, AbstractString}
    record.refID < 0 && return nothing
    return @inbounds record.refs.names[record.refID + 1]
end

function next_reference_name(record::BAMRecord)::Union{Nothing, AbstractString}
    record.next_refID < 0 && return nothing
    return @inbounds record.refs.names[record.next_refID + 1]
end

"""
    references(record::BAMRecord)::References

Get the reference sequences that the references of `record` refer to.
Records read from a `BAMReader` share the reader's references.

See also: [`reference_index`](@ref), [`next_reference_index`](@ref)
"""
references(record::BAMRecord) = record.refs

"""
    reference_index(record::BAMRecord)::Union{Nothing, Int}

Return the 1-based index of the reference of `record` in [`references(record)`](@ref references),
or `nothing` if the record has no reference.

This is the reference ID (refID) of the BAM format plus one.

# Examples
```julia
index = reference_index(record)
if index !== nothing
    @assert references(record)[index].name == reference_name(record)
end
```

See also: [`reference`](@ref XAM2.reference), [`reference_name`](@ref), [`next_reference_index`](@ref)
"""
function reference_index(record::BAMRecord)::Union{Nothing, Int}
    record.refID < 0 && return nothing
    return (record.refID % Int) + 1
end

"""
    next_reference_index(record::BAMRecord)::Union{Nothing, Int}

Return the 1-based index of the reference of the next segment in the template (RNEXT)
in [`references(record)`](@ref references), or `nothing` if it has no reference.

This is the next reference ID (next_refID) of the BAM format plus one.

See also: [`reference_index`](@ref), [`next_reference_name`](@ref)
"""
function next_reference_index(record::BAMRecord)::Union{Nothing, Int}
    record.next_refID < 0 && return nothing
    return (record.next_refID % Int) + 1
end

"""
    XAM2.reference(record::BAMRecord)::Union{Nothing, Reference}

Return the reference of `record`, with its name and length, or `nothing` if
the record has no reference.

See also: [`reference_index`](@ref), [`XAM2.Reference`](@ref)
"""
function reference(record::BAMRecord)::Union{Nothing, Reference}
    record.refID < 0 && return nothing
    return @inbounds record.refs[record.refID + 1]
end

function reference_length(record::BAMRecord)::Union{Int, Nothing}
    record.refID < 0 && return nothing
    return @inbounds(record.refs.lengths[record.refID + 1]) % Int
end

# Invariant: The variable-length fields fit in the data.
# Offsets are computed as Int, so they can't overflow.
cigar_offset(record::BAMRecord) = record.l_read_name % Int
seq_offset(record::BAMRecord) = cigar_offset(record) + 4 * (record.n_cigar_op % Int)
qual_offset(record::BAMRecord) = seq_offset(record) + cld(record.l_seq % Int, 2)
aux_offset(record::BAMRecord) = qual_offset(record) + record.l_seq % Int

# The 1-based position of the first `=` in the sequence, or `nothing` if there is none
function first_eq_seq_pos(record::BAMRecord)::Union{Nothing, Int}
    l_seq = record.l_seq % Int
    span = (seq_offset(record) + 1):qual_offset(record)
    mem = @inbounds ImmutableMemoryView(record.data)[span]
    # Each byte holds two bases, the first in the high nibble. `=` is encoded as 0.
    # If the length is odd, the low nibble of the last byte is padding.
    n_whole_bytes = l_seq >> 1
    # Fast path: Skip chunks of bytes with no `=`
    i = 1
    chunk_size = 128
    @inbounds while i <= n_whole_bytes - chunk_size + 1
        any = false
        for j in i:(i + chunk_size - 1)
            any |= iszero(mem[j] & 0x0f) | iszero(mem[j] & 0xf0)
        end
        any && break
        i += chunk_size
    end
    @inbounds while i <= n_whole_bytes
        byte = mem[i]
        iszero(byte & 0xf0) && return 2 * i - 1
        iszero(byte & 0x0f) && return 2 * i
        i += 1
    end
    if isodd(l_seq) && iszero(@inbounds(mem[end]) & 0xf0)
        return l_seq
    end
    return nothing
end

"""
    sequence(rec::BAMRecord)::Union{Nothing, BAMSeqString}

Get the query sequence of `rec` as a [`BAMSeqString`](@ref), or
`nothing` if the sequence is not available.

The string is constructed in constant time, and is a view into the record's data, so it
is not separately allocated. It is never empty.

With BioSequences loaded, use `XAM2.try_sequence(BioSequence, rec)` for a mutable
view as a `BioSequence`, or `try_sequence(LongDNA{4}, rec)` for a copy.

!!! warning
    The returned sequence aliases the record. Mutating the record will
    invalidate the sequence
"""
function sequence(record::BAMRecord)::Union{BAMSeqString, Nothing}
    iszero(record.l_seq) && return nothing
    span = (seq_offset(record) + 1):qual_offset(record)
    return BAMSeqString(@inbounds(ImmutableMemoryView(record.data)[span]), record.l_seq % Int)
end

"""
    quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}

Return the PHRED quality scores of `record`, or `nothing` if the quality is
missing (stored as all `0xff`, or the sequence is missing). The returned scores are
never empty.

The scores have no offset, so e.g. a PHRED score of 31 is `UInt8(31)`, as opposed to
`UInt8(31 + 33)` in SAM. They are stored like this in BAM, so this is a view of the
record's data.

To modify the scores, use [`set_quality_scores!`](@ref), or
[`XAM2.unsafe_quality_bytes`](@ref).

See also: [`quality_string`](@ref)

!!! warning
    The returned view aliases the record. Mutating the record will
    invalidate the view.
"""
function quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    qual = unsafe_quality_bytes(record)
    # Missing quality is encoded as all-0xFF
    # A missing sequence has no quality bytes, and its quality is also missing
    (isempty(qual) || @inbounds(first(qual)) == 0xff) && return nothing
    return ImmutableMemoryView(qual)
end

"""
    quality_string(record::BAMRecord)::Union{Nothing, StringView{QualityView}}

Return the quality of `record` as a PHRED+33-encoded string, as in the QUAL field of
SAM, or `nothing` if the quality is missing (see [`quality_scores`](@ref)).
The returned string is never empty.

The string is backed by a lazy [`XAM2.QualityView`](@ref), which adds 33 to the
PHRED scores of the record on access. Use `String` to materialize it.

!!! warning
    The returned string aliases the record. Mutating the record will
    invalidate the string.
"""
function quality_string(record::BAMRecord)::Union{Nothing, StringView{QualityView}}
    qual = quality_scores(record)
    return qual === nothing ? nothing : StringView(bam_to_sam_quality(qual))
end

"""
    aux_data(record::BAMRecord)::Union{XAMAuxData.BAM.MutableAuxiliary, XAMError}

Return the auxiliary data as a mutable `XAMAuxData.BAM.Auxiliary`, or an `XAMError`
of kind `XAMErrors.invalid_aux_data` if it is malformed, i.e. if the tags or the
boundaries of the fields can't be found.
If the record does not contain any auxiliary fields, the returned object is empty.

This checks the structure of all fields, which is O(n) in the size of the auxiliary data.
The values are not validated: Loading an invalid value (e.g. a non-numeric integer) returns
an `XAMAuxData.Error`. Check all values with `isvalid`.
To skip the check, see [`aux_data_unchecked`](@ref).

!!! warning
    The returned `Auxiliary` object aliases the record.
    Mutating the record will invalidate the auxiliary data.

See the documentation of `XAMAuxData.BAM.Auxiliary` for more details.
"""
function aux_data(record::BAMRecord)::Union{XAMAuxData.BAM.MutableAuxiliary, XAMError}
    aux = aux_data_unchecked(record)
    return is_well_formed(aux) ? aux : XAMError(XAMErrors.invalid_aux_data)
end

"""
    aux_data_unchecked(record::BAMRecord)::XAMAuxData.BAM.MutableAuxiliary

Like [`aux_data`](@ref), but without checking that the auxiliary data is well-formed.
This avoids scanning all fields, e.g. when only a few fields are accessed.
If the auxiliary data is malformed, iterating it or loading values may throw an
exception, or return wrong values. Check with `XAMAuxData.is_well_formed`.

!!! warning
    The returned `Auxiliary` object aliases the record.
    Mutating the record will invalidate the auxiliary data.
"""
function aux_data_unchecked(record::BAMRecord)::XAMAuxData.BAM.MutableAuxiliary
    return XAMAuxData.BAM.Auxiliary(record.data, aux_offset(record) + 1)
end

const CIGAR_LUT = let
    u = UInt64(0)
    for (i, byte) in enumerate(b"MIDNSHP=X")
        u |= (byte % UInt64) << (7 * i - 7)
    end
    u
end

# Invariant: pos and next_pos are in BAM_POSITIONS
function left_position(record::BAMRecord)::Union{Int, Nothing}
    return record.pos < 0 ? nothing : (record.pos % Int) + 1
end

function next_position(record::BAMRecord)::Union{Int, Nothing}
    return record.next_pos < 0 ? nothing : (record.next_pos % Int) + 1
end

function mapping_quality(record::BAMRecord)::Union{Int, Nothing}
    return record.mapq == 0xff ? nothing : record.mapq % Int
end

function template_length(record::BAMRecord)::Union{Int, Nothing}
    return iszero(record.tlen) ? nothing : record.tlen % Int
end

# Errors are thrown by outlined functions, since constructing an error inline makes
# the caller set up a GC frame. See the TODO about this in XAM2.jl.
@noinline throw_invalid_value(msg::String, value) = throw(ArgumentError(string(msg, ", got ", repr(value))))
@noinline throw_argument_error(msg::String) = throw(ArgumentError(msg))
@noinline throw_key_error(key) = throw(KeyError(key))

# Setters for the fixed-length fields of `BAMRecord`.
# Each setter accepts exactly the values its getter can return, plus `nothing` for
# missing, so that `getter(setter!(record, x)) == x`. Invalid input throws, and leaves
# the record unchanged.

"""
    set_flags!(record::BAMRecord, flags::SAMFlags)::BAMRecord

Set the flags of `record`, and return `record`.
Throw an `ArgumentError` if `flags` lacks the `UNMAP` flag, and the record is not
placed (see [`set_placement!`](@ref)), since a mapped record must have a reference and
position. To map an unplaced record, set its placement first.

See also: [`flags`](@ref)
"""
function set_flags!(record::BAMRecord, flags::SAMFlags)
    if !(SAMFlag.UNMAP in flags) && !is_placed(record)
        throw_argument_error("Cannot clear the UNMAP flag of a record without a reference and position")
    end
    record.flags = flags
    return record
end

is_placed(record::BAMRecord) = (record.refID >= 0) & (record.pos >= 0)

"""
    set_mapping_quality!(record::BAMRecord, mapq::Union{Nothing, Integer})::BAMRecord

Set the mapping quality (MAPQ) of `record` to `mapq`, which must be in `0:254`,
or to missing if `mapq` is `nothing`, and return `record`.
Throw an `ArgumentError` if `mapq` is out of range. In particular, to set a missing
mapping quality, pass `nothing`, not 255.

See also: [`mapping_quality`](@ref)
"""
function set_mapping_quality!(record::BAMRecord, mapq::Union{Nothing, Integer})
    record.mapq = if mapq === nothing
        0xff
    else
        in(mapq, 0:254) || throw_invalid_value("Mapping quality must be in 0:254", mapq)
        mapq % UInt8
    end
    return record
end

"""
    set_placement!(record::BAMRecord, ref::Union{Integer, AbstractString}, pos::Integer)::BAMRecord
    set_placement!(record::BAMRecord, nothing)::BAMRecord

Set the reference (RNAME) and 1-based leftmost position (POS) of `record`, and return
`record`. With `nothing`, remove both, making the record unplaced.

`ref` may be a 1-based index into [`references(record)`](@ref references), like the one
returned by [`reference_index`](@ref), or the name of a reference.
`pos` must be within the reference. The end of the alignment is not checked against
the length of the reference, since it depends on the CIGAR.

Throw, leaving the record unchanged:
* A `BoundsError` if `ref` is an out-of-bounds index, or a `KeyError` if no
  reference has the name `ref`.
* An `ArgumentError` if `pos` is not in `1:len`, where `len` is the length of the
  reference.
* An `ArgumentError` when removing the placement of a record that does not have the
  `UNMAP` flag, since a mapped record must have a reference and position, or of a
  record with a CIGAR of more than 65535 operations, which BAM can only store on a
  record with a reference.

See also: [`reference_index`](@ref), [`left_position`](@ref), [`set_next_placement!`](@ref)
"""
function set_placement!(record::BAMRecord, ref::Union{Integer, AbstractString}, pos::Integer)
    ref_id = bam_reference_id(record.refs, ref)
    record.pos = bam_position(record.refs, ref_id, pos)
    record.refID = ref_id
    return record
end

function set_placement!(record::BAMRecord, ::Nothing)
    if !(SAMFlag.UNMAP in record.flags)
        throw_argument_error("Cannot remove the reference and position of a mapped record")
    end
    if has_placeholder_cigar(record)
        throw_argument_error("Cannot remove the reference of a record with more than 65535 CIGAR operations")
    end
    record.refID = Int32(-1)
    record.pos = Int32(-1)
    return record
end

# If the CIGAR field is the placeholder for a CIGAR stored in the CG tag
function has_placeholder_cigar(record::BAMRecord)::Bool
    span = (cigar_offset(record) + 1):seq_offset(record)
    return is_placeholder_cigar(record, @inbounds ImmutableMemoryView(record.data)[span])
end

"""
    set_next_placement!(record::BAMRecord, ref::Union{Integer, AbstractString}, pos::Integer)::BAMRecord
    set_next_placement!(record::BAMRecord, nothing)::BAMRecord

Set the reference (RNEXT) and 1-based leftmost position (PNEXT) of the next segment in
the template of `record`, and return `record`. With `nothing`, remove both.
`ref` and `pos` are as in [`set_placement!`](@ref), and throw likewise.

The flags are not checked or changed. Make sure that when setting or clearing placement,
the MUNMAP and MREVERSE flags are correctly set.

See also: [`next_reference_index`](@ref), [`next_position`](@ref)
"""
function set_next_placement!(record::BAMRecord, ref::Union{Integer, AbstractString}, pos::Integer)
    ref_id = bam_reference_id(record.refs, ref)
    record.next_pos = bam_position(record.refs, ref_id, pos)
    record.next_refID = ref_id
    return record
end

function set_next_placement!(record::BAMRecord, ::Nothing)
    record.next_refID = Int32(-1)
    record.next_pos = Int32(-1)
    return record
end

# The 0-based BAM position of a 1-based position.
# Invariant: The result is in `BAM_POSITIONS`.
# `pos` must be within the reference with the 0-based ID `ref_id`.
# Invariant: The result is in `BAM_POSITIONS`, since reference lengths are at most
# `typemax(Int32)`.
function bam_position(refs::References, ref_id::Int32, pos::Integer)::Int32
    ref_len = @inbounds refs.lengths[ref_id + 1]
    in(pos, 1:ref_len) || throw_position_out_of_range(pos, ref_len)
    return Int32(pos) - Int32(1)
end

@noinline function throw_position_out_of_range(pos::Integer, ref_len::Int32)
    throw(ArgumentError(lazy"Position must be in 1:$ref_len, the length of the reference, got $pos"))
end

# The 0-based reference ID of `ref` in `refs`.
# Invariant: The result is in `0:(length(refs) - 1)`.
function bam_reference_id(refs::References, ref::Union{Integer, AbstractString})::Int32
    ref isa AbstractString && return @something reference_id(refs, ref) throw_key_error(ref)
    in(ref, 1:length(refs)) || throw(BoundsError(refs, ref))
    # References has at most typemax(Int32) elements
    return Int32(ref) - Int32(1)
end

"""
    set_template_length!(record::BAMRecord, tlen::Union{Nothing, Integer})::BAMRecord

Set the observed template length (TLEN) of `record` to `tlen`, or to missing if `tlen`
is `nothing`, and return `record`.
Throw an `ArgumentError` if `tlen` is zero, or does not fit in an `Int32`.
To set a missing template length, pass `nothing`, not zero.

See also: [`template_length`](@ref)
"""
function set_template_length!(record::BAMRecord, tlen::Union{Nothing, Integer})
    record.tlen = if tlen === nothing
        Int32(0)
    else
        if iszero(tlen) || !in(tlen, typemin(Int32):typemax(Int32))
            throw_invalid_value("Template length must be a nonzero Int32", tlen)
        end
        Int32(tlen)
    end
    return record
end

# Aliasing
# Setters that change the length of a variable-length field move the bytes after it.
# A setter's input may be a view of the record's own data, e.g. when setting the name
# from an auxiliary string value with `set_query_name!(r, aux_data(r)["RX"])`, in which
# case the input's bytes may be moved or overwritten before they are read.
# For views of predictable types (see `DataView`), setters handle this without copying:
# 1. `data_indices` finds the indices in `data` of the input's bytes, if it is a view of it.
#    String inputs are checked by their `codeunits`, which for `StringView`s of
#    `MemoryView`s, and substrings of those, is a `MemoryView`. This covers the strings
#    returned by the accessors, so setters only need to handle arrays.
# 2. `replace_region!` writes the new content at a point where the input is intact, and
#    passes the input's current location to the writing function.
# Aliasing by other types, e.g. user-defined views, is not detected, and is the user's
# responsibility, so the docstrings of such setters must warn about it.

# SubArrays with fast linear indexing are contiguous when their first index is an
# `AbstractUnitRange`, or all their indices are scalars.
const ContiguousByteSubArray{N, P} = SubArray{UInt8, N, P, I, true} where {
    I <: Union{Tuple{AbstractUnitRange, Vararg{Any}}, Tuple{Vararg{Integer}}},
}

# Types of byte vectors that are checked for aliasing the record's data.
# Views returned by accessors of `BAMRecord`, and the code units of strings returned
# by them, are `MemoryView`s.
const DataView = Union{
    MemoryView{UInt8},
    ContiguousByteSubArray{1, <:Union{Vector{UInt8}, Memory{UInt8}, MemoryView{UInt8}}},
}

# The indices in `data` of the bytes of `x`, if `x` is a non-empty view of the memory of
# `data`, else `nothing`. This only compares addresses, so `x` may be a view of any memory.
data_indices(::Vector{UInt8}, ::Any) = nothing

function data_indices(data::Vector{UInt8}, x::DataView)::Union{Nothing, UnitRange{Int}}
    n = length(x) % UInt
    iszero(n) && return nothing
    # Comparing addresses does not dereference them, so the objects need not be preserved
    p = UInt(pointer(x))
    d = UInt(pointer(data))
    p < d && return nothing
    start = p - d
    len = length(data) % UInt
    (start <= len) & (n <= len - start) || return nothing
    return ((start + 1) % Int):((start + n) % Int)
end

# Replace the `old_len` bytes of `record.data` after the 0-based `offset` with `new_len`
# bytes written by `write!(data, src)`, moving the bytes after the region.
# Throw an `ArgumentError` if the record would be too long for BAM (invariant 1),
# in which case the record is unchanged and `write!` is not called. This is checked
# before modifying anything.
#
# `src` is the input of the setter, as given by `data_indices` (or `nothing` if it
# does not alias `data`). `write!` gets the index in `data` where the input's bytes
# currently start, or `nothing`. It must write `new_len` bytes at `offset + 1`, and it
# must copy aliased input with `unsafe_move_bytes!`, since the input may overlap the
# region. `write!` is called when the input is intact:
# * When shrinking, before moving the trailing bytes, so the input has not moved.
# * When growing, after moving them, since there is no room before. Input before or
#   inside the region is not moved, and input after it is moved by the growth.
#   Input that overlaps the end of the region would be split in two by the move;
#   this is an internal error, and setters must rule it out beforehand.
#
# This breaks the invariants of the record: The caller must update the length fields
# (e.g. `l_read_name`) to match, and `write!` must not throw.
@inline function replace_region!(
        write!::F,
        record::BAMRecord,
        offset::Int,
        old_len::Int,
        new_len::Int,
        src::Union{Nothing, UnitRange{Int}},
    )::Nothing where {F}
    data = record.data
    old_total = length(data)
    # TODO: Remove this assert when the package is ready for release.
    # Its throw is outlined, since an inline `@assert` makes the function set up a GC frame.
    if !((offset >= 0) & (old_len >= 0) & (new_len >= 0) & (offset + old_len <= old_total))
        throw_assertion_error("invalid region in replace_region!")
    end
    new_total = old_total - old_len + new_len
    new_total > typemax(UInt32) - BAM_FIXED_BYTES && throw_record_too_long()
    n_trailing = old_total - offset - old_len
    region_end = offset + old_len
    trailing_src = region_end + 1
    trailing_dst = offset + new_len + 1
    if new_len <= old_len
        write!(data, src === nothing ? nothing : first(src))
        if new_len < old_len
            unsafe_move_bytes!(data, trailing_dst, trailing_src, n_trailing)
            resize!(data, new_total)
        end
    else
        moved_src = if src === nothing
            nothing
        elseif first(src) > region_end
            first(src) + (new_len - old_len)
        else
            # TODO: Remove this assert when the package is ready for release.
            last(src) <= region_end || throw_assertion_error("input split by replace_region!")
            first(src)
        end
        resize!(data, new_total)
        unsafe_move_bytes!(data, trailing_dst, trailing_src, n_trailing)
        write!(data, moved_src)
    end
    return nothing
end

# Copy `n` bytes of `data` from index `src` to index `dst`, without bounds checking.
# The regions may overlap, since the copy is a `memmove`.
# This uses pointers, because slicing `MemoryView`s adds checks for empty slices,
# which the compiler can't remove.
@inline function unsafe_move_bytes!(data::Vector{UInt8}, dst::Int, src::Int, n::Int)::Nothing
    GC.@preserve data unsafe_copyto!(pointer(data, dst), pointer(data, src), n % UInt)
    return nothing
end

@noinline throw_assertion_error(msg::String) = throw(AssertionError(msg))

@noinline function throw_record_too_long()
    throw_argument_error("Record would be too long to be represented in BAM")
end

"""
    set_query_name!(record::BAMRecord, name::Union{Nothing, AbstractString})::BAMRecord

Set the query name (QNAME) of `record` to `name`, or to missing if `name` is `nothing`,
and return `record`.
Throw an `ArgumentError` if `name` is not a valid name, i.e. if it does not match the
regex `[!-?A-~]{1,254}` from the SAM specification, or is `"*"`, which denotes
a missing name. To set a missing name, pass `nothing`.

This moves the rest of the record's data, which is O(n) in the size of the record.

!!! warning
    `name` must not alias the record's data, except for strings returned by this package,
    such as [`query_name(record)`](@ref query_name), a string value from
    [`aux_data(record)`](@ref aux_data), and substrings of these, which are handled.

See also: [`query_name`](@ref)
"""
function set_query_name!(record::BAMRecord, name::Union{Nothing, AbstractString})
    if name === nothing
        unsafe_set_query_name!(record, "*")
    else
        utf8 = to_utf8_string(name)
        is_valid_query_name(utf8) || throw_invalid_value("Invalid query name", name)
        unsafe_set_query_name!(record, utf8)
    end
    return record
end

# String types whose code units are UTF-8 bytes. For other string types, code units
# may be wider than a byte (e.g. UTF-16), or a byte in the ASCII range may not be that
# ASCII character (e.g. UTF-7), so they must be converted before handling their bytes.
const UTF8String = Union{String, SubString{String}, StringView, SubString{<:StringView}}

# Return `s` as a `UTF8String`: `s` itself if it is one, else a converted `String`.
# The conversion is outlined, since it allocates, and other string types are rare.
to_utf8_string(s::UTF8String) = s
@noinline to_utf8_string(s::AbstractString)::String = String(s)

# Set the query name to `name`, which must be a valid name, or "*".
# Strings are passed to this and `is_valid_query_name` instead of their code units,
# since passing a `CodeUnits` to a function that is not inlined makes the caller set up
# a GC frame.
function unsafe_set_query_name!(record::BAMRecord, name::UTF8String)::Nothing
    bytes = codeunits(name)
    n = length(bytes)
    # A view of the data can't overlap the end of the name region, which is the name's
    # null byte, since a valid name has no null bytes. So `replace_region!` can handle it.
    src = data_indices(record.data, bytes)
    replace_region!(record, 0, record.l_read_name % Int, n + 1, src) do data, src_start
        if src_start === nothing
            for i in 1:n
                @inbounds data[i] = bytes[i]
            end
        else
            unsafe_move_bytes!(data, 1, src_start, n)
        end
        # Invariant 3: The name is null terminated
        @inbounds data[n + 1] = 0x00
    end
    record.l_read_name = (n + 1) % UInt8
    return nothing
end

# Check if `name` matches `[!-?A-~]{1,254}`, and is not `*`.
# A valid name therefore has no null bytes, and its length plus the null byte fits in
# `l_read_name`.
function is_valid_query_name(name::UTF8String)::Bool
    bytes = codeunits(name)
    in(length(bytes), 1:254) || return false
    length(bytes) == 1 && @inbounds(bytes[1]) == UInt8('*') && return false
    for c in bytes
        (in(c, 0x21:0x7e) & (c != 0x40)) || return false
    end
    return true
end

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, XAMError}
    encoding = @something get_cigar_encoding(record) return nothing
    cig = CIGARStrings.try_parse(BAMCIGAR, encoding)
    cig isa CIGARError && return XAMError(cig)
    return check_query_length(cig, record.l_seq % Int)
end

# CIGARs with more than 65535 operations are stored in the CG auxiliary tag,
# with a placeholder CIGAR `kSmN`, where k is the sequence length, and m the
# reference length of the alignment.
@inline function is_placeholder_cigar(record::BAMRecord, cigar_encoding::ImmutableMemoryView{UInt8})
    length(cigar_encoding) == 8 || return false
    record.refID < 0 && return false
    u1 = unsafe_load_integer_le(UInt32, cigar_encoding, 1)
    u2 = unsafe_load_integer_le(UInt32, cigar_encoding, 5)
    # Like htslib, don't check m, which is the reference length of the alignment,
    # and so is only known after loading the CIGAR.
    return ((u1 & 0x0f) == 0x04) & ((u2 & 0x0f) == 0x03) & ((u1 >> 4) == record.l_seq)
end

function get_cigar_encoding(x::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    iszero(x.n_cigar_op) && return nothing
    span = (cigar_offset(x) + 1):seq_offset(x)
    vw = @inbounds ImmutableMemoryView(x.data)[span]
    return is_placeholder_cigar(x, vw) ? load_aux_cigar(x, vw) : vw
end

# Get the CIGAR in the CG tag, or return the placeholder CIGAR if the tag is
# not present, like htslib. Outlined, because it's expected to be rare.
@noinline function load_aux_cigar(
        x::BAMRecord,
        placeholder::ImmutableMemoryView{UInt8},
    )::ImmutableMemoryView{UInt8}
    aux = aux_data(x)
    aux isa XAMError && return placeholder
    v = get(aux, "CG", nothing)
    return v isa AbstractVector{UInt32} ? parent(v) : placeholder
end

"""
    set_cigar!(record::BAMRecord, cigar::Union{Nothing, BAMCIGAR})::BAMRecord

Set the CIGAR of `record` to `cigar`, or to missing if `cigar` is `nothing`, and return
`record`. To set a `CIGARStrings.CIGAR`, convert it with `BAMCIGAR(cigar)` first.

The CIGAR must match the sequence: If `record` has a sequence, the query length of
`cigar` must equal the length of the sequence. To change both, use
[`set_cigar_and_sequence!`](@ref).

!!! warning
    This removes any `CG` auxiliary field of `record`. A `CG` field is only meaningful
    as the storage of the record's CIGAR, so it would otherwise be stale.

A CIGAR with more than 65535 operations does not fit in BAM's CIGAR field, so it is
stored in a new `CG` auxiliary field, with a placeholder in the CIGAR field, as described
in the SAM specification. [`cigar`](@ref) returns the full CIGAR.

Throw an `ArgumentError`, leaving the record unchanged, if:
* `cigar` has no operations. To set a missing CIGAR, pass `nothing`.
* `record` has a sequence, and the query length of `cigar` differs from its length.
* The record would be too long to be represented in BAM.
* `cigar` has more than 65535 operations, and the record has no reference, or the
  length of the sequence or the reference length of `cigar` is at least `2^28`, or the
  record's auxiliary data is malformed. In these cases the CIGAR can't be stored.

This moves the rest of the record's data, which is O(n) in the size of the record.
`cigar` may be a view of the record's own data, e.g. `cigar(record)`.

See also: [`cigar`](@ref), [`set_cigar_and_sequence!`](@ref)
"""
function set_cigar!(record::BAMRecord, cig::Union{Nothing, BAMCIGAR})
    if cig === nothing
        old_len = 4 * (record.n_cigar_op % Int)
        replace_region!((_, _) -> nothing, record, cigar_offset(record), old_len, 0, nothing)
        record.n_cigar_op = 0x0000
    else
        n_ops = length(cig)
        iszero(n_ops) && throw_argument_error("Cannot set an empty CIGAR. To set a missing CIGAR, pass `nothing`")
        check_cigar_query_length(cig, record.l_seq % Int)
        n_ops > typemax(UInt16) && return set_long_cigar!(record, cig)
        unsafe_set_cigar_field!(record, MemoryView(cig))
    end
    delete_cg_tag!(record)
    return record
end

function check_cigar_query_length(cig::BAMCIGAR, l_seq::Int)::Nothing
    query_len = CIGARStrings.query_length(cig)
    if !iszero(l_seq) && query_len != l_seq
        throw_query_length_mismatch(query_len, l_seq)
    end
    return nothing
end

@noinline function throw_query_length_mismatch(query_len::Int, l_seq::Int)
    throw(ArgumentError(lazy"Query length of CIGAR, $query_len, differs from sequence length $l_seq"))
end

# Replace the CIGAR field with `bytes`, the encoding of at most 65535 CIGAR operations.
# `bytes` may alias the record's data. The CG tag is not changed.
function unsafe_set_cigar_field!(record::BAMRecord, bytes::ImmutableMemoryView{UInt8})::Nothing
    offset = cigar_offset(record)
    old_len = 4 * (record.n_cigar_op % Int)
    new_len = length(bytes)
    src = data_indices(record.data, bytes)
    # Aliased input that overlaps the end of the region would be split by the move
    # when growing. This is not possible for input from accessors, so it's outlined.
    region_end = offset + old_len
    if src !== nothing && (new_len > old_len) & (first(src) <= region_end) & (last(src) > region_end)
        return unsafe_set_cigar_field!(record, copy(bytes))
    end
    replace_region!(record, offset, old_len, new_len, src) do data, src_start
        @inline
        if src_start === nothing
            GC.@preserve data bytes unsafe_copyto!(pointer(data, offset + 1), pointer(bytes), new_len % UInt)
        else
            unsafe_move_bytes!(data, offset + 1, src_start, new_len)
        end
    end
    record.n_cigar_op = (new_len >>> 2) % UInt16
    return nothing
end


# Delete the CG tag, which is stale after a new CIGAR is set.
# If the auxiliary data is malformed, the tag can't be found, but it doesn't need to
# be deleted either: `load_aux_cigar` never reads malformed auxiliary data.
function delete_cg_tag!(record::BAMRecord)::Nothing
    aux = aux_data(record)
    aux isa XAMError || delete!(aux, "CG")
    return nothing
end

# Set a CIGAR with more than 65535 operations, which is stored in the CG tag, with the
# placeholder `kSmN` in the CIGAR field (see `store_cigar_placeholder!`).
# Outlined, since this is rare. The query length of `cig` has already been checked.
@noinline function set_long_cigar!(record::BAMRecord, cig::BAMCIGAR)::BAMRecord
    l_seq = record.l_seq % Int
    check_long_cigar(record, cig, l_seq)
    old_len = 4 * (record.n_cigar_op % Int)
    check_long_cigar_fits(record, cig, length(record.data) - old_len + 8)
    # The edits below move the record's data, so aliased input is copied, e.g. if the
    # input is the record's current CIGAR stored in the CG tag
    bytes = unaliased(record.data, MemoryView(cig))
    offset = cigar_offset(record)
    ref_len = CIGARStrings.ref_length(cig)
    replace_region!(record, offset, old_len, 8, nothing) do data, _
        store_cigar_placeholder!(data, offset, l_seq, ref_len)
    end
    record.n_cigar_op = 0x0002
    replace_cg_tag!(record, bytes)
    return record
end

# Throw if a CIGAR with more than 65535 operations can't be stored on `record`, if
# the record's sequence length is `l_seq`.
function check_long_cigar(record::BAMRecord, cig::BAMCIGAR, l_seq::Int)::Nothing
    # The placeholder is only recognised on records with a reference
    if record.refID < 0
        throw_argument_error("A CIGAR with more than 65535 operations can only be set on a record with a reference")
    end
    if max(l_seq, CIGARStrings.ref_length(cig)) > 0x0fff_ffff
        throw_argument_error("The lengths of a CIGAR with more than 65535 operations must be less than 2^28")
    end
    # A CG tag in malformed auxiliary data would not be read by `load_aux_cigar`
    if aux_data(record) isa XAMError
        throw_argument_error("A CIGAR with more than 65535 operations can't be set on a record with malformed auxiliary data")
    end
    return nothing
end

# Throw if the record would be too long after appending a CG tag holding `cig`, where
# `len` is the length of the data before appending. Checked before modifying anything.
# This ignores the removal of an existing CG tag, so it may reject records that would
# barely fit.
function check_long_cigar_fits(record::BAMRecord, cig::BAMCIGAR, len::Int)::Nothing
    new_total = len + CG_HEADER_BYTES + 4 * length(cig)
    new_total > typemax(UInt32) - BAM_FIXED_BYTES && throw_record_too_long()
    return nothing
end

# Store the placeholder CIGAR `kSmN` after the 0-based `offset`, where k is the sequence
# length, and m the reference length (see `is_placeholder_cigar`).
function store_cigar_placeholder!(data::Vector{UInt8}, offset::Int, l_seq::Int, ref_len::Int)::Nothing
    mem = MemoryView(data)
    inbounds_store_le!(mem, offset + 1, ((l_seq % UInt32) << 4) | 0x04) # S
    inbounds_store_le!(mem, offset + 5, ((ref_len % UInt32) << 4) | 0x03) # N
    return nothing
end

# Replace any CG tag with one holding `bytes`, the encoding of a CIGAR.
# The auxiliary data must be well-formed, and `bytes` must not alias the record's data.
function replace_cg_tag!(record::BAMRecord, bytes::ImmutableMemoryView{UInt8})::Nothing
    delete!(aux_data_unchecked(record), "CG")
    copyto!(append_cg_header!(record.data, length(bytes) >>> 2), bytes)
    return nothing
end

# Return `x`, or a copy of it if it aliases `data`.
# The copy is outlined, since it allocates, and aliasing is rare.
@inline function unaliased(data::Vector{UInt8}, x::ImmutableMemoryView{UInt8})::ImmutableMemoryView{UInt8}
    return data_indices(data, x) === nothing ? x : copy_outlined(x)
end

@noinline copy_outlined(x) = copy(x)

# Bytes of a CG tag before its array: tag, type 'B', element type 'I', and count
const CG_HEADER_BYTES = 8

# Append the header of a CG:B:I tag with `n_ops` elements to `data`, grow `data` to fit
# the elements, and return a view of the bytes of the elements, to be filled by the caller.
function append_cg_header!(data::Vector{UInt8}, n_ops::Int)::MutableMemoryView{UInt8}
    old_len = length(data)
    resize!(data, old_len + CG_HEADER_BYTES + 4 * n_ops)
    mem = MemoryView(data)[(old_len + 1):end]
    copyto!(mem, b"CGBI")
    inbounds_store_le!(mem, 5, n_ops % UInt32)
    return mem[(CG_HEADER_BYTES + 1):end]
end

# Sequence and quality setters.
# The sequence input is normalized to `nothing`, a `BAMSeqString` whose packed bytes are
# copied, or a `UTF8String`, whose bytes are encoded. Their code units are a
# `MemoryView` for `StringView`s and their substrings, so they are checked for aliasing.
# Strings are passed instead of their code units, since passing a `CodeUnits` to a
# function that is not inlined allocates (see `unsafe_set_query_name!`).
# The quality input is normalized to `nothing`, an `ImmutableMemoryView` if possible
# (so it can be checked for aliasing, and copied without bounds checks), or else left as is.
normalize_seq(::Nothing) = nothing
normalize_seq(seq::BAMSeqString) = seq
normalize_seq(seq::AbstractString) = to_utf8_string(seq)

normalize_qual(::Nothing) = nothing
normalize_qual(qual::Union{DataView, Vector{UInt8}, Memory{UInt8}}) = ImmutableMemoryView(qual)
normalize_qual(qual::AbstractVector{UInt8}) = qual

seq_length(::Nothing) = 0
seq_length(seq::BAMSeqString) = ncodeunits(seq)
seq_length(seq::UTF8String) = ncodeunits(seq)

# Check the sequence and quality input, and return the sequence length.
function check_seq_qual(seq, qual::Union{Nothing, AbstractVector{UInt8}})::Int
    if seq === nothing
        qual === nothing || throw_argument_error("A record without a sequence can't have quality scores")
        return 0
    end
    n = seq_length(seq)
    iszero(n) && throw_argument_error("Cannot set an empty sequence. To set a missing sequence, pass `nothing`")
    if seq isa UTF8String
        for byte in codeunits(seq)
            # Branchless: The mask keeps the index in bounds, and the comparison
            # rejects the bytes it changes
            ok = (byte < 0x80) & (@inbounds(STRICT_ENCODE_LUT[(byte & 0x7f) + 1]) != 0xff)
            ok || throw_invalid_value("Sequence can only contain the characters `=ACMGRSVTWYHKDBN` and their lowercase", Char(byte))
        end
    end
    qual === nothing || check_qual(qual, n)
    return n
end

# Check that `qual` are `n` quality scores which SAM can represent, i.e. in `0:93`,
# the printable characters `!` to `~`.
function check_qual(qual::AbstractVector{UInt8}, n::Int)
    length(qual) == n || throw_qual_length_mismatch(length(qual), n)
    bad = false
    for q in qual
        bad |= q > 0x5d
    end
    bad && throw_bad_qual(qual)
    return nothing
end

@noinline function throw_bad_qual(qual::AbstractVector{UInt8})
    for q in qual
        q > 0x5d && throw_invalid_value("Quality scores must be in 0:93", Int(q))
    end
    return unreachable()
end

@noinline function throw_qual_length_mismatch(qual_len::Int, n::Int)
    throw(ArgumentError(lazy"Length of quality scores, $qual_len, differs from sequence length $n"))
end

# Return the inputs, copying those that alias `data`.
function unaliased_seq_qual(data::Vector{UInt8}, seq, qual)
    seq = if seq isa BAMSeqString
        data_indices(data, seq.data) === nothing ? seq : BAMSeqString(copy_outlined(seq.data), seq.len)
    elseif seq isa UTF8String && data_indices(data, codeunits(seq)) !== nothing
        StringView(copy_outlined(ImmutableMemoryView(codeunits(seq))))
    else
        seq
    end
    qual = qual isa ImmutableMemoryView{UInt8} ? unaliased(data, qual) : qual
    return (seq, qual)
end

# Number of bytes of the encoded sequence and quality of length `n`
seq_qual_bytes(n::Int) = cld(n, 2) + n

# Write the sequence of length `n` and its quality after the 0-based `offset`.
# The input must have been checked, and must not alias `data`.
function write_seq_qual!(data::Vector{UInt8}, offset::Int, seq, qual, n::Int)::Nothing
    seq === nothing && return nothing
    mem = MemoryView(data)
    n_seq_bytes = cld(n, 2)
    seq_mem = @inbounds mem[(offset + 1):(offset + n_seq_bytes)]
    if seq isa BAMSeqString
        # Invariant of `BAMSeqString`: Its data has `n_seq_bytes` bytes
        @inbounds copyto!(seq_mem, seq.data)
    else
        # Like `bamseq!`, but the input has already been checked
        fill!(seq_mem, 0x00)
        bytes = codeunits(seq)
        for i in 1:n
            byte = @inbounds bytes[i]
            encoding = @inbounds STRICT_ENCODE_LUT[byte + 1]
            @inbounds seq_mem[(i + 1) >>> 1] |= encoding << ((4 * isodd(i)) & 7)
        end
    end
    qual_mem = @inbounds mem[(offset + n_seq_bytes + 1):(offset + n_seq_bytes + n)]
    if qual === nothing
        # Missing quality is encoded as all 0xff
        fill!(qual_mem, 0xff)
    else
        # The length of `qual` was checked to be `n`
        @inbounds copyto!(qual_mem, qual)
    end
    return nothing
end

"""
    set_sequence!(
        record::BAMRecord,
        seq::Union{Nothing, AbstractString},
        qual::Union{Nothing, AbstractVector{UInt8}},
    )::BAMRecord

Set the sequence (SEQ) of `record` to `seq`, and its quality scores (QUAL) to `qual`,
and return `record`. `seq` may be `nothing` for a missing sequence, and `qual` may be
`nothing` for missing quality scores.

`seq` may contain the characters `=ACMGRSVTWYHKDBN`, i.e. `=` and the IUPAC nucleotide
codes, which are the characters BAM can store. Lowercase letters are accepted, and
stored as uppercase.
`qual` holds PHRED scores without offset, like the return value of
[`quality_scores`](@ref).

If `record` has a CIGAR, the length of the sequence can't be changed, since it must
match the CIGAR. To change both, use [`set_cigar_and_sequence!`](@ref).

Throw an `ArgumentError`, leaving the record unchanged, if:
* `seq` is empty, or contains other characters than `=ACMGRSVTWYHKDBN` and their
  lowercase. To set a missing sequence, pass `nothing`.
* `qual` is not `nothing`, and `seq` is `nothing`, or `qual` has a different length
  than `seq`, or `qual` contains a score above 93, which SAM can't represent.
* `record` has a CIGAR, and the length of `seq` differs from the current sequence length
  (including setting a sequence on a record with a missing sequence, or vice versa).
* The record would be too long to be represented in BAM (about 2 GB).

This moves the rest of the record's data, if the length changes.

See also: [`sequence`](@ref), [`quality_scores`](@ref), [`set_cigar_and_sequence!`](@ref)
"""
function set_sequence!(
        record::BAMRecord,
        seq::Union{Nothing, AbstractString},
        qual::Union{Nothing, AbstractVector{UInt8}},
    )
    s = normalize_seq(seq)
    q = normalize_qual(qual)
    n = check_seq_qual(s, q)
    l_seq = record.l_seq % Int
    if !iszero(record.n_cigar_op) && n != l_seq
        throw_argument_error("Cannot change the sequence length of a record with a CIGAR. Use `set_cigar_and_sequence!`")
    end
    # New names, since reassigned captured variables are boxed
    (seq_in, qual_in) = unaliased_seq_qual(record.data, s, q)
    offset = seq_offset(record)
    replace_region!(record, offset, seq_qual_bytes(l_seq), seq_qual_bytes(n), nothing) do data, _
        write_seq_qual!(data, offset, seq_in, qual_in, n)
    end
    record.l_seq = n % UInt32
    return record
end

"""
    set_cigar_and_sequence!(
        record::BAMRecord,
        cigar::Union{Nothing, BAMCIGAR},
        seq::Union{Nothing, AbstractString},
        qual::Union{Nothing, AbstractVector{UInt8}},
    )::BAMRecord

Set the CIGAR, sequence and quality scores of `record` together, and return `record`.
This is needed to change the length of the sequence of a record with a CIGAR, e.g.
when trimming a read. The arguments are as in [`set_cigar!`](@ref) and
[`set_sequence!`](@ref), and if both `cigar` and `seq` are given, the query length of
`cigar` must equal the length of `seq`.

Throws like [`set_cigar!`](@ref) and [`set_sequence!`](@ref), except that the sequence
length may change. This removes any `CG` auxiliary field, like `set_cigar!`.

This moves the rest of the record's data once, which is O(n) in the size of the record.

See also: [`set_cigar!`](@ref), [`set_sequence!`](@ref)
"""
function set_cigar_and_sequence!(
        record::BAMRecord,
        cig::Union{Nothing, BAMCIGAR},
        seq::Union{Nothing, AbstractString},
        qual::Union{Nothing, AbstractVector{UInt8}},
    )
    s = normalize_seq(seq)
    q = normalize_qual(qual)
    n = check_seq_qual(s, q)
    (seq_in, qual_in) = unaliased_seq_qual(record.data, s, q)
    return replace_cigar_and_seq_qual!(record, cig, n) do data, offset
        write_seq_qual!(data, offset, seq_in, qual_in, n)
    end
end

# Replace the CIGAR, sequence and quality of `record` with `cig`, and `n` bases written
# by `write_seq_qual!(data, offset)` after the 0-based `offset`, moving the rest of the
# data once. `cig` is checked against `n` before anything is modified.
# `write_seq_qual!` is called before the CIGAR is written, when the old CIGAR, sequence
# and quality are still intact at their old offsets (see `replace_region!`), so it may
# read them, as long as it writes nothing before the old sequence. It must not throw.
function replace_cigar_and_seq_qual!(
        write_seq_qual!::F,
        record::BAMRecord,
        cig::Union{Nothing, BAMCIGAR},
        n::Int,
    )::BAMRecord where {F}
    data = record.data
    old_len = 4 * (record.n_cigar_op % Int) + seq_qual_bytes(record.l_seq % Int)
    # Variables captured by the closure below are assigned only once, since
    # reassigned captured variables are boxed
    is_long = cig !== nothing && length(cig) > typemax(UInt16)
    if cig !== nothing
        iszero(length(cig)) && throw_argument_error("Cannot set an empty CIGAR. To set a missing CIGAR, pass `nothing`")
        check_cigar_query_length(cig, n)
        if is_long
            check_long_cigar(record, cig, n)
            check_long_cigar_fits(record, cig, length(data) - old_len + 8 + seq_qual_bytes(n))
        end
    end
    # Copied if aliased, since `write_seq_qual!` may overwrite the old CIGAR
    cig_bytes = cig === nothing ? nothing : unaliased(data, MemoryView(cig))
    ref_len = is_long ? CIGARStrings.ref_length(something(cig)) : 0
    cig_len = cig_bytes === nothing ? 0 : is_long ? 8 : length(cig_bytes)
    offset = cigar_offset(record)
    replace_region!(record, offset, old_len, cig_len + seq_qual_bytes(n), nothing) do data, _
        write_seq_qual!(data, offset + cig_len)
        if is_long
            store_cigar_placeholder!(data, offset, n, ref_len)
        elseif cig_bytes !== nothing
            GC.@preserve data cig_bytes unsafe_copyto!(pointer(data, offset + 1), pointer(cig_bytes), cig_len % UInt)
        end
    end
    record.n_cigar_op = (cig_len >>> 2) % UInt16
    record.l_seq = n % UInt32
    if is_long
        replace_cg_tag!(record, something(cig_bytes))
    else
        delete_cg_tag!(record)
    end
    return record
end

"""
    slice_sequence!(record::BAMRecord, idx::AbstractUnitRange{<:Integer})::BAMRecord
    slice_sequence!(record::BAMRecord, idx::AbstractUnitRange{<:Integer}, cigar::Union{Nothing, BAMCIGAR})::BAMRecord

Keep only the bases at the 1-based indices `idx` of the sequence of `record`, and their
quality scores, and return `record`. This is equivalent to
`set_sequence!(record, sequence(record)[idx], quality_scores(record)[idx])`, but works
in place, without allocating. Missing quality scores stay missing.

With a CIGAR argument, the CIGAR is set like with [`set_cigar_and_sequence!`](@ref),
which is needed to change the sequence length of a record with a CIGAR. E.g. to hard clip
a read with CIGAR `2S8M`, use `slice_sequence!(record, 3:10, BAMCIGAR(CIGAR("2H8M")))`.
A missing CIGAR (`nothing`) is valid for any record, also mapped ones. This
function doesn't change the position, so if the alignment start moves, use
[`set_placement!`](@ref).

Throw an `ArgumentError`, leaving the record unchanged, if:
* `idx` is empty, or not in `1:n` for a sequence of length `n`. A record with a missing
  sequence has no bases to keep. To set a missing sequence, use [`set_sequence!`](@ref)
  or [`set_cigar_and_sequence!`](@ref).
* No CIGAR is given, `record` has a CIGAR, and `idx` is not the whole sequence.
* A CIGAR is given, and [`set_cigar!`](@ref) would throw with a sequence of length
  `length(idx)`.

Like `set_cigar_and_sequence!`, the three-argument form removes any `CG` auxiliary field.

This moves the rest of the record's data, which is O(n) in the size of the record.

See also: [`set_sequence!`](@ref), [`set_cigar_and_sequence!`](@ref)
"""
function slice_sequence!(record::BAMRecord, idx::AbstractUnitRange{<:Integer})
    (start, n) = check_slice(record, idx)
    if !iszero(record.n_cigar_op) && n != record.l_seq
        throw_argument_error("Cannot change the sequence length of a record with a CIGAR. Pass the new CIGAR to `slice_sequence!`")
    end
    l_seq = record.l_seq % Int
    offset = seq_offset(record)
    qual_missing = is_qual_missing(record)
    replace_region!(record, offset, seq_qual_bytes(l_seq), seq_qual_bytes(n), nothing) do data, _
        unsafe_slice_seq_qual!(data, offset, l_seq, start, n, qual_missing)
    end
    record.l_seq = n % UInt32
    return record
end

function slice_sequence!(record::BAMRecord, idx::AbstractUnitRange{<:Integer}, cig::Union{Nothing, BAMCIGAR})
    (start, n) = check_slice(record, idx)
    l_seq = record.l_seq % Int
    old_offset = seq_offset(record)
    qual_missing = is_qual_missing(record)
    return replace_cigar_and_seq_qual!(record, cig, n) do data, offset
        unsafe_slice_seq_qual!(data, old_offset, l_seq, start, n, qual_missing)
        # Move the sliced bases and quality from the old to the new sequence offset
        unsafe_move_bytes!(data, offset + 1, old_offset + 1, seq_qual_bytes(n))
    end
end

# Check that `idx` is a nonempty range of indices of the sequence of `record`, and
# return its start and length.
function check_slice(record::BAMRecord, idx::AbstractUnitRange{<:Integer})::Tuple{Int, Int}
    start = Int(first(idx))::Int
    stop = Int(last(idx))::Int
    stop < start && throw_argument_error("Cannot slice to an empty sequence")
    if (start < 1) | (stop > record.l_seq)
        throw(ArgumentError(lazy"Range $idx is out of bounds for sequence of length $(record.l_seq)"))
    end
    return (start, stop - start + 1)
end

is_qual_missing(record::BAMRecord) = quality_scores(record) === nothing

# Keep the `n` bases starting at the 1-based index `start` of the sequence of length
# `l_seq` and their quality, moving them to the start of the sequence after the 0-based
# `offset`. Their quality follows them, and is set to all `0xff` if `qual_missing`.
# Writes nothing outside the old sequence and quality, which must be intact, and
# `start:(start + n - 1)` must be in `1:l_seq`.
@inline function unsafe_slice_seq_qual!(
        data::Vector{UInt8},
        offset::Int,
        l_seq::Int,
        start::Int,
        n::Int,
        qual_missing::Bool,
    )::Nothing
    n_bytes = cld(n, 2)
    # Base `i` is in byte `(i + 1) >> 1`, in the high nibble if `i` is odd.
    # New byte `j` is read from old bytes at index `j` or higher, so the bytes can be
    # moved in order, without overwriting bytes that are yet to be read.
    if isodd(start)
        # The bases are byte-aligned, so the bytes are moved as they are
        unsafe_move_bytes!(data, offset + 1, offset + (start + 1) >> 1, n_bytes)
    else
        # New byte `j` is the low nibble of old byte `k`, and the high nibble of byte
        # `k + 1`, if the latter is in the sequence.
        k0 = offset + start >> 1
        n_full = n >> 1
        for j in 1:n_full
            k = k0 + j - 1
            @inbounds data[offset + j] = (data[k] << 4) | (data[k + 1] >> 4)
        end
        if isodd(n)
            @inbounds data[offset + n_bytes] = data[k0 + n_full] << 4
        end
    end
    # Zero the padding nibble, which may hold a base of the old sequence
    if isodd(n)
        @inbounds data[offset + n_bytes] &= 0xf0
    end
    qual_dst = offset + n_bytes + 1
    if qual_missing
        @inbounds fill!(view(data, qual_dst:(qual_dst + n - 1)), 0xff)
    else
        unsafe_move_bytes!(data, qual_dst, offset + cld(l_seq, 2) + start, n)
    end
    return nothing
end

# In-place operations on the sequence and quality scores, which change no lengths

"""
    unsafe_sequence_bytes(record::BAMRecord)::MutableMemoryView{UInt8}

Return a mutable view of the encoded sequence of `record`, which is empty if the
sequence is missing. This is a low-level operation for modifying bases in place.
For a mutable biological sequence, use `try_sequence(BioSequence, record)`, which needs
BioSequences to be loaded.

The view has `cld(n, 2)` bytes, for a sequence of length `n`. Each byte holds two bases
as 4-bit codes, the first in the high nibble. The code of a base is its 0-based index in
`"=ACMGRSVTWYHKDBN"`. If `n` is odd, the low nibble of the last byte is padding, which
must be zero.

This is unsafe in that the view is unchecked: It's the caller's responsibility to keep
the padding nibble zero. Any change to the record that moves its sequence, e.g. setting
its name or CIGAR, invalidates the view, after which writing to it corrupts the record.

# Examples
```julia
# Set the 3rd base, which is in the high nibble of the 2nd byte, to `T`
bytes = XAM2.unsafe_sequence_bytes(record)
bytes[2] = (bytes[2] & 0x0f) | (0x08 << 4)
```

See also: [`sequence`](@ref), [`XAM2.unsafe_quality_bytes`](@ref)
"""
function unsafe_sequence_bytes(record::BAMRecord)::MutableMemoryView{UInt8}
    # If the sequence is missing, this span is empty
    span = (seq_offset(record) + 1):qual_offset(record)
    return @inbounds MemoryView(record.data)[span]
end

"""
    unsafe_quality_bytes(record::BAMRecord)::MutableMemoryView{UInt8}

Return a mutable view of the encoded quality of `record`. This is a low-level operation
for modifying quality scores in place. For a checked alternative, use
[`set_quality_scores!`](@ref).

The view has one byte per base, holding the PHRED score without offset, so it is empty
if the sequence is missing. Unlike [`quality_scores`](@ref), this returns the bytes even
if the quality is missing, in which case they are all `0xff`.

This is unsafe in that the view is unchecked: Writing a score above 93 makes the record
impossible to write as SAM, and writing `0xff` to the first byte makes the quality
read as missing. Any change to the record that moves its quality, e.g. setting its name,
CIGAR or sequence, invalidates the view, after which writing to it corrupts the record.

# Examples
```julia
# Bin quality scores in place
bytes = XAM2.unsafe_quality_bytes(record)
if !isempty(bytes) && first(bytes) != 0xff
    map!(q -> q < 20 ? 0x0a : 0x1e, bytes, bytes)
end
```

See also: [`XAM2.unsafe_sequence_bytes`](@ref)
"""
function unsafe_quality_bytes(record::BAMRecord)::MutableMemoryView{UInt8}
    # If the sequence is missing, this span is empty
    span = (qual_offset(record) + 1):aux_offset(record)
    return @inbounds MemoryView(record.data)[span]
end

"""
    set_quality_scores!(record::BAMRecord, qual::Union{Nothing, AbstractVector{UInt8}})::BAMRecord

Set the quality scores of `record` to `qual`, or mark them as missing if `qual` is
`nothing`. Like [`quality_scores`](@ref), `qual` holds PHRED scores without offset.
This works whether or not the quality is currently missing, and does not move any data
of the record, since the length of the quality is the length of the sequence.
To change the length, use [`set_sequence!`](@ref).

Throw an `ArgumentError`, leaving the record unchanged, if:
* `record` has no sequence
* `qual` has a different length than the sequence, or contains a score above 93,
  which SAM can't represent.

`qual` may alias the record's quality, e.g. be its own [`quality_scores`](@ref).

# Examples
```julia
# Bin quality scores, reusing a buffer across records
buffer = UInt8[]
for record in records
    scores = quality_scores(record)
    scores === nothing && continue
    resize!(buffer, length(scores))
    map!(q -> q < 20 ? 0x0a : 0x1e, buffer, scores)
    set_quality_scores!(record, buffer)
end
```

See also: [`XAM2.unsafe_quality_bytes`](@ref)
"""
function set_quality_scores!(record::BAMRecord, qual::Union{Nothing, AbstractVector{UInt8}})
    bytes = unsafe_quality_bytes(record)
    isempty(bytes) && throw_argument_error("A record without a sequence can't have quality scores")
    if qual === nothing
        fill!(bytes, 0xff)
    else
        check_qual(qual, length(bytes))
        copyto!(bytes, qual)
    end
    return record
end

"""
    set_record!(
        record::BAMRecord;
        query_name::Union{Nothing, AbstractString} = nothing,
        flags::SAMFlags = SAMFlags([SAMFlag.UNMAP]),
        reference::Union{Nothing, Integer, AbstractString} = nothing,
        left_position::Union{Nothing, Integer} = nothing,
        mapping_quality::Union{Nothing, Integer} = nothing,
        cigar::Union{Nothing, BAMCIGAR} = nothing,
        sequence::Union{Nothing, AbstractString} = nothing,
        quality_scores::Union{Nothing, AbstractVector{UInt8}} = nothing,
        next_reference::Union{Nothing, Integer, AbstractString} = nothing,
        next_position::Union{Nothing, Integer} = nothing,
        template_length::Union{Nothing, Integer} = nothing,
    )::BAMRecord

Set all fields of `record` except the auxiliary data, and return `record`.
This is the way to build a record from scratch, e.g. in a read simulator.
To avoid allocating a new record for each read, create one record with
[`BAMRecord(refs)`](@ref BAMRecord), and call `set_record!` on it for each read.

The fields that are not given are set to those of an empty record, so the result does
not depend on the previous content of `record`: Missing values, no CIGAR or sequence,
and flags with only `UNMAP` set. The auxiliary data is removed, and can be set
afterwards with [`aux_data`](@ref).

The keywords are named after the getters, and accept the same values as the
corresponding setters, e.g. [`set_placement!`](@ref) for `reference` and
`left_position`, and [`set_cigar_and_sequence!`](@ref) for `cigar`, `sequence` and
`quality_scores`. `reference` and `left_position` must be given together, and so must
`next_reference` and `next_position`.

Unlike when calling the individual setters, the order of the fields doesn't matter:
E.g. `flags` may lack `UNMAP` if `reference` and `left_position` are given.

Throws like the individual setters. If this throws, `record` is valid, but its content
is unspecified. Inputs may alias the record's data, e.g. be its own
[`query_name`](@ref), in which case they are copied.

# Examples
```julia
record = BAMRecord(refs)
for read in simulated_reads
    set_record!(
        record;
        query_name = read.name,
        flags = SAMFlags([SAMFlag.REVERSE]),
        reference = "chr1",
        left_position = read.pos,
        mapping_quality = 60,
        cigar = read.cigar,
        sequence = read.seq,
        quality_scores = read.qual,
    )
    aux_data(record)["NM"] = read.nm
    write(writer, record)
end
```

See also: [`reset!`](@ref)
"""
function set_record!(
        record::BAMRecord;
        query_name::Union{Nothing, AbstractString} = nothing,
        flags::SAMFlags = SAMFlags(0x0004),
        reference::Union{Nothing, Integer, AbstractString} = nothing,
        left_position::Union{Nothing, Integer} = nothing,
        mapping_quality::Union{Nothing, Integer} = nothing,
        cigar::Union{Nothing, BAMCIGAR} = nothing,
        sequence::Union{Nothing, AbstractString} = nothing,
        quality_scores::Union{Nothing, AbstractVector{UInt8}} = nothing,
        next_reference::Union{Nothing, Integer, AbstractString} = nothing,
        next_position::Union{Nothing, Integer} = nothing,
        template_length::Union{Nothing, Integer} = nothing,
    )
    check_both_or_neither(reference, left_position, "reference", "left_position")
    check_both_or_neither(next_reference, next_position, "next_reference", "next_position")
    # Copy inputs that alias the data, since `reset!` overwrites it.
    # Strings are normalized first, so their code units are checked for aliasing.
    data = record.data
    name = unaliased_string(data, query_name === nothing ? nothing : to_utf8_string(query_name))
    cig = if cigar === nothing || data_indices(data, MemoryView(cigar)) === nothing
        cigar
    else
        copy_outlined(cigar)
    end
    (seq, qual) = unaliased_seq_qual(data, normalize_seq(sequence), normalize_qual(quality_scores))
    reset!(record)
    # Placement must be set before flags without `UNMAP`, and before a long CIGAR
    set_query_name!(record, name)
    reference === nothing || set_placement!(record, reference, something(left_position))
    set_flags!(record, flags)
    set_mapping_quality!(record, mapping_quality)
    next_reference === nothing || set_next_placement!(record, next_reference, something(next_position))
    set_template_length!(record, template_length)
    set_cigar_and_sequence!(record, cig, seq, qual)
    return record
end

function check_both_or_neither(a, b, name_a::String, name_b::String)::Nothing
    if (a === nothing) != (b === nothing)
        throw(ArgumentError(lazy"`$name_a` and `$name_b` must be given together"))
    end
    return nothing
end

# Return `s`, or a copy of it if its code units alias `data`.
unaliased_string(::Vector{UInt8}, ::Nothing) = nothing

function unaliased_string(data::Vector{UInt8}, s::UTF8String)
    data_indices(data, codeunits(s)) === nothing && return s
    return StringView(copy_outlined(ImmutableMemoryView(codeunits(s))))
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
    XAM2.try_read!(reader::BAMReader)::Union{BAMRecord, Nothing, XAMError}

Attempt to read the next record of `reader`, similar to iterating the reader.
If the reader is EOF, return `nothing`.
If the next record is malformed, return an `XAMError` with the virtual offset
of the start of the record.
Once this function returns an `XAMError`, subsequent calls to this function
will return the same error until the reader is closed or seeked.

The structure of the record is validated, e.g. that its fields fit in the record,
and that its reference IDs and positions are in bounds. The content of variable-length
fields, e.g. the CIGAR, is validated when it is accessed.

This function does not throw if the reader's data is malformed.
"""
function try_read!(x::BAMReader)::Union{BAMRecord, Nothing, XAMError}
    return try_read_record!(x, false)
end

# Like `try_read!`, but if `new_record` is set, read into a new record instead of the
# reader's record. This copies the data once, whereas copying the reader's record
# would copy it twice.
function try_read_record!(x::BAMReader, new_record::Bool)::Union{BAMRecord, Nothing, XAMError}
    isnothing(x.error) || return x.error

    # Check if EOF and return that error if so
    buffer = try_get_nonempty_buffer(x.reader)
    # The malformed block is where the next record starts, so seeking there reproduces
    # the error.
    buffer isa BGZFError && return fail_read!(x, XAMError(buffer; offset = block_start(buffer)))
    isnothing(buffer) && return nothing
    # Since the buffer is non-empty, this is in the block containing the record's start
    offset = virtual_position(x.reader)

    record = x.record

    # Read the block size, and the fixed-length fields into the record
    memview = fill_bytes_to_data!(x, 4 + BAM_FIXED_BYTES)
    memview isa ImmutableMemoryView || return fail_read!(x, fill_error(memview, offset))
    stored_block_size = unsafe_load_integer_le(UInt32, memview, 1)
    unsafe_load_fixed_fields(record, @inbounds memview[5:end])
    # Clear undefined flag bits, which `SAMFlags` must not have
    record.flags = SAMFlags(Integer(record.flags))
    err = check_fixed_fields(record, stored_block_size)
    err === nothing || return fail_read!(x, XAMError(err; offset))

    data_len = (stored_block_size - UInt32(BAM_FIXED_BYTES)) % Int
    if new_record
        # Allocate after checking the fixed fields, so obviously malformed records are
        # rejected without allocating. A malformed block size may still pass the check
        # and allocate up to 4 GiB, but the pages are only touched as data is read.
        new = new_uninitialized_record(x.refs, Vector{UInt8}(undef, data_len))
        copy_fixed_fields!(new, record)
        # The reader's record may have been returned by `try_read!`, so it must
        # satisfy the structural invariants of `BAMRecord`.
        reset!(record)
        record = new
    end

    # Read all the variable length fields. For the reader's record, this overwrites the
    # fixed field bytes, which we have already moved into the integer fields of the struct
    memview = fill_bytes_to_data!(x.reader, record.data, data_len)
    memview isa ImmutableMemoryView || return fail_read!(x, fill_error(memview, offset))
    resize!(record.data, data_len)
    @inbounds(record.data[record.l_read_name]) == 0x00 ||
        return fail_read!(x, XAMError(XAMErrors.invalid_read_name; offset))
    return record
end

# Set the reader's sticky error, and reset its record, since a failed read
# may leave the record in a state that violates its invariants.
@noinline function fail_read!(x::BAMReader, err::XAMError)::XAMError
    reset!(x.record)
    x.error = err
    return err
end

# The error when a record could not be read in full, where `offset` is the start of
# the record.
function fill_error(result::Union{Nothing, BGZFError}, offset::VirtualOffset)::XAMError
    return if result === nothing
        XAMError(XAMErrors.truncated_file; offset)
    else
        XAMError(result; offset)
    end
end

# The virtual offset of the start of the block of a BGZF error, if known
function block_start(e::BGZFError)::Union{Nothing, VirtualOffset}
    file_offset = e.file_offset
    return file_offset === nothing ? nothing : VirtualOffset(file_offset, 0)
end

# Check the structural invariants of `record` that only depend on its fixed fields,
# see the comment in the definition of `BAMRecord`.
# Valid 0-based positions in BAM, where -1 means missing
const BAM_POSITIONS = Int32(-1):(typemax(Int32) - Int32(1))

function check_fixed_fields(record::BAMRecord, stored_block_size::UInt32)::Union{Nothing, XAMErrorKind}
    # The block size includes the fixed fields after the block size field.
    stored_block_size < BAM_FIXED_BYTES && return XAMErrors.invalid_block_length
    data_len = (stored_block_size - UInt32(BAM_FIXED_BYTES)) % Int
    iszero(record.l_read_name) && return XAMErrors.invalid_read_name
    n_refs = length(record.refs)
    is_valid_reference_id(record.refID, n_refs) || return XAMErrors.invalid_reference_id
    is_valid_reference_id(record.next_refID, n_refs) || return XAMErrors.invalid_reference_id
    in(record.pos, BAM_POSITIONS) || return XAMErrors.invalid_position
    in(record.next_pos, BAM_POSITIONS) || return XAMErrors.invalid_next_position
    # Computed as Int, so it can't overflow
    l_seq = record.l_seq % Int
    min_len = record.l_read_name % Int + 4 * (record.n_cigar_op % Int) + cld(l_seq, 2) + l_seq
    min_len > data_len && return XAMErrors.invalid_block_length
    return nothing
end

"""
    mutable struct BAMWriter{T <: AbstractBufWriter}

Writer of BAM files.

    BAMWriter(
        io::Union{IO, AbstractBufWriter},
        header::SAMHeader,
        refs::References = References(header);
        n_workers::Int = min(4, Threads.nthreads()),
        compress_level::Int = 6,
    )::BAMWriter

    BAMWriter(bgzf::BGZFWriter, header::SAMHeader, refs::References = References(header))

Construct a `BAMWriter` writing to `io`, by wrapping `io` in a `BGZFLib.BGZFWriter`.
The header and references are written immediately.
Closing the resulting writer will close `io`.
To obtain the `BGZFWriter` without closing it, use `BufferIO.take_inner(::BAMWriter)`.

The text of `header` is written verbatim. If `header` contains `@SQ` lines, they must
match `refs`, else an `ArgumentError` is thrown. If `refs` is not passed, it is
created from the `@SQ` lines of `header`.

When writing a record, its reference IDs are translated from the record's references
to the writer's references by reference name. Translation is skipped if the record
and the writer share the same `References` object, e.g. when constructing the writer
with `BAMWriter(io, get_header(reader), references(reader))`, or if the two
`References` have identical content.
A record cannot be written if its reference is not present in the writer's
references, or if the reference has a different length in the writer.
See [`XAM2.try_write`](@ref).

# Examples
```julia
BAMReader(open("in.bam")) do reader
    # The text header of a BAM file may be malformed, see `get_header`
    header = get_header(reader)
    header isa XAMError && throw(header)
    BAMWriter(open("out.bam", "w"), header, references(reader)) do writer
        for record in reader
            is_mapped(record) && write(writer, record)
        end
    end
end
```
"""
mutable struct BAMWriter{T <: AbstractBufWriter}
    const writer::BGZFWriter{T}
    const refs::References
    # A Record stores a Reference, but what do we do if the Reference of the record
    # does not match the one in this writer?
    # Solution: Translate the reference, then cache the translated reference in this
    # field. This amortizes translation of writes of many reads from a different Reference.
    # A map of `nothing` means the source has the same content as `refs`.
    # The `last_*` fields is a single-element cache in front of `maps`.
    last_src::Union{Nothing, References}
    last_map::Union{Nothing, Vector{Int32}}
    const maps::IdDict{References, Union{Nothing, Vector{Int32}}}
    # Set by `close` and `take_inner`. After `take_inner`, `writer` is still open,
    # but no longer owned by this writer.
    closed::Bool
end

references(writer::BAMWriter) = writer.refs

BAMRecord(io::Union{BAMReader, BAMWriter}) = BAMRecord(references(io))

function BAMWriter(
        io::Union{AbstractBufWriter, IO},
        header::SAMHeader,
        refs::References = References(header);
        n_workers::Int = min(4, Threads.nthreads()),
        compress_level::Int = 6,
    )
    # Check before creating the BGZF writer, so we don't leak it if this throws
    check_writer_header(header, refs)
    bgzf = BGZFWriter(io; n_workers, compress_level)
    return BAMWriter(bgzf, header, refs)
end

function BAMWriter(
        writer::BGZFWriter{T},
        header::SAMHeader,
        refs::References = References(header),
    ) where {T}
    check_writer_header(header, refs)
    header_text = MemoryView(header)

    # Write magic bytes BAM\1
    write_le(writer, 0x014d4142)
    write_le(writer, length(header_text) % UInt32)
    write(writer, header_text)
    write_le(writer, length(refs) % UInt32)

    # Write reference sequences
    for ref in refs
        cu = codeunits(ref.name)
        write_le(writer, (length(cu) + 1) % UInt32)
        write(writer, cu)
        write(writer, 0x00)
        write_le(writer, ref.length % UInt32)
    end

    maps = IdDict{References, Union{Nothing, Vector{Int32}}}()
    return BAMWriter{T}(writer, refs, nothing, nothing, maps, false)
end

# Write `x` as 4 little-endian bytes.
function write_le(writer::AbstractBufWriter, x::UInt32)::Nothing
    buffer = get_nonempty_buffer(writer, 4)
    buffer === nothing && return write_le_slowpath(writer, x)
    inbounds_store_le!(buffer, 1, x)
    @inbounds consume(writer, 4)
    return nothing
end

# The writer can't buffer 4 bytes at once, so write them one at a time
@noinline function write_le_slowpath(writer::AbstractBufWriter, x::UInt32)::Nothing
    for shift in 0:8:24
        write(writer, (x >>> shift) % UInt8)
    end
    return nothing
end

function check_writer_header(header::SAMHeader, refs::References)::Nothing
    if length(MemoryView(header)) > typemax(Int32)
        throw(ArgumentError("Header must be 2147483647 bytes or shorter"))
    end
    check_header_references(header, refs)
    return nothing
end

function BAMWriter(f::Function, args::Vararg{Any, N}; kwargs...) where {N}
    writer = BAMWriter(args...; kwargs...)
    return try
        f(writer)
    finally
        close(writer)
    end
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

# Closing twice is a no-op, so that closing after `take_inner` does not close the
# taken writer, e.g. in the `do`-block form of `BAMWriter`.
function Base.close(writer::BAMWriter)
    writer.closed && return nothing
    writer.closed = true
    return close(writer.writer)
end

function Base.flush(writer::BAMWriter)
    writer.closed && throw(IOError(IOErrorKinds.ClosedIO))
    return flush(writer.writer)
end

Base.isopen(writer::BAMWriter) = !writer.closed && isopen(writer.writer)

"""
    BufferIO.take_inner(writer::BAMWriter{T})::BGZFWriter{T}

Close `writer` and take its `BGZFLib.BGZFWriter`, without closing or flushing it.
This is the `BGZFWriter` `writer` was constructed from, or the one it created when
constructed from an `io`.
Records are written directly to the `BGZFWriter`, so it contains everything written
to `writer`. Closing the `BGZFWriter` completes the BAM file by writing the BGZF EOF block,
unless it was created with `append_empty = false`. It can also be passed to
`BufferIO.take_inner` to get the underlying writer without closing it.

Calling this on an already closed `writer` throws a `BufferIO.IOError` of kind `ClosedIO`.
After this, writing to `writer` throws the same error, and `close(writer)` does nothing.

# Examples
```julia
buffer = IOBuffer()
writer = BAMWriter(buffer, header)
write(writer, record)
bgzf = BufferIO.take_inner(writer)
# Write the remaining data and the EOF block, without closing `buffer`
flush(BufferIO.take_inner(bgzf))
bam_bytes = take!(buffer)
```
"""
function BufferIO.take_inner(writer::BAMWriter)
    writer.closed && throw(IOError(IOErrorKinds.ClosedIO))
    writer.closed = true
    return writer.writer
end

"""
    XAM2.try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMError}

Write `record` to `writer`, and return the number of bytes written.
If the record can't be written, write nothing, and return an `XAMError`. This happens if:
* The reference or next reference of `record` is not present in `writer`,
  or has a different length in `writer`.
* The CIGAR of `record` is malformed, or doesn't match the sequence length,
  see [`cigar`](@ref).
Errors from the underlying IO are thrown.

The bin field, which is used for indexing, is computed from the position, CIGAR and
flags of `record` when writing, so it is correct even if the record has been edited.
The auxiliary data is not validated.

See also: [`BAMWriter`](@ref)
"""
function try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMError}
    result = try_write_kind(writer, record)
    result isa XAMErrorKind && return XAMError(result)
    return result
end

# The BAI bin of a record with 0-based position `pos` (-1 if missing), whose CIGAR spans
# `ref_len` reference positions (0 if it has no CIGAR). Following the SAM specification,
# section 5.3, and htslib, unmapped records, and records whose CIGAR consumes no reference
# positions, are binned as spanning one position, and records with no position get
# `reg2bin(-1, 0)`, which is 4680.
function compute_bin(pos::Int32, ref_len::Int, flags::SAMFlags)::UInt16
    pos < 0 && return UInt16(4680)
    len = SAMFlag.UNMAP in flags ? 1 : max(ref_len, 1)
    # As `UInt64`, since the end may not fit in 32 bits
    beg = pos % UInt64
    return reg2bin(beg, beg + len % UInt64)
end

function try_write_kind(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMError, XAMErrorKind}
    # The stored `bin` is not trusted, since the record may have been edited.
    # Computing it needs the CIGAR, which is validated in the same pass.
    cig = cigar(record)
    cig isa XAMError && return cig
    ref_len = cig === nothing ? 0 : CIGARStrings.ref_length(cig)
    bin = compute_bin(record.pos, ref_len, record.flags)
    src = record.refs
    map = src === writer.refs ? nothing : reference_map!(writer, src)
    (ref_id, next_ref_id) = if map === nothing
        n_refs = length(writer.refs)
        is_valid_reference_id(record.refID, n_refs) || return XAMErrors.invalid_reference_id
        is_valid_reference_id(record.next_refID, n_refs) || return XAMErrors.invalid_reference_id
        (record.refID, record.next_refID)
    else
        ref_id = translate_reference_id(map, record.refID)
        ref_id isa XAMErrorKind && return ref_id
        next_ref_id = translate_reference_id(map, record.next_refID)
        next_ref_id isa XAMErrorKind && return next_ref_id
        (ref_id, next_ref_id)
    end
    return write_record(writer, record, ref_id, next_ref_id, bin)
end

function Base.write(writer::BAMWriter, record::BAMRecord)
    result = try_write(writer, record)
    result isa XAMError && throw(result)
    return result
end

is_valid_reference_id(id::Int32, n_refs::Int) = (id ≥ -1) & (id < n_refs)

function reference_map!(writer::BAMWriter, src::References)::Union{Nothing, Vector{Int32}}
    src === writer.last_src && return writer.last_map
    map = get!(() -> build_reference_map(src, writer.refs), writer.maps, src)
    writer.last_src = src
    writer.last_map = map
    return map
end

# Write the record, but with the given reference IDs and bin instead of the record's own.
function write_record(
        writer::BAMWriter,
        record::BAMRecord,
        ref_id::Int32,
        next_ref_id::Int32,
        bin::UInt16,
    )::Int
    writer.closed && throw(IOError(IOErrorKinds.ClosedIO))
    bgzf = writer.writer
    buffer = get_buffer(bgzf)
    n_header_bytes = 4 + BAM_FIXED_BYTES
    if length(buffer) < n_header_bytes
        grow_buffer(bgzf)
        buffer = get_buffer(bgzf)
        # This is unreachable, because calling `grow_buffer` on a BGZFWriter
        # will always make sure the next buffer is a full BGZF block,
        # which is 64 KiB.
        length(buffer) < n_header_bytes && error("unreachable")
    end
    inbounds_store_le!(buffer, 1, block_size(record))
    GC.@preserve record buffer begin
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            (pointer(buffer) + 4)::Ptr{Cchar},
            src::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    # Byte offsets of refID, bin and next_refID, after the block size
    inbounds_store_le!(buffer, 5, ref_id)
    @inbounds buffer[15] = bin % UInt8
    @inbounds buffer[16] = (bin >> 8) % UInt8
    inbounds_store_le!(buffer, 25, next_ref_id)
    @inbounds consume(bgzf, n_header_bytes)
    # Write variable bytes
    return n_header_bytes + write(bgzf, record.data)
end

@inline function inbounds_store_le!(mem::MutableMemoryView{UInt8}, i::Int, x::Union{Int32, UInt32})
    u = x % UInt32
    for offset in 0:3
        @inbounds mem[i + offset] = (u >> (8 * offset)) % UInt8
    end
    return nothing
end

# The ID of the reference named `name` in `refs`, or -1 if `name` is `nothing`,
# i.e. the record has no reference.
function sam_reference_id(refs::References, name::Union{Nothing, AbstractString})::Union{Int32, XAMErrorKind}
    name === nothing && return Int32(-1)
    return @something reference_id(refs, name) XAMErrors.unknown_reference
end

"""
    XAM2.try_copy!(dst::BAMRecord, src::SAMRecord)::Union{BAMRecord, XAMError}
    XAM2.try_copy!(dst::SAMRecord, src::BAMRecord)::Union{SAMRecord, XAMError}

Convert `src` to the format of `dst`, overwriting `dst`, and return `dst`.
`copy!(dst, src)` does the same, but throws the error.

If `src` can't be converted, `dst` is reset to an empty record, and an `XAMError`
is returned. This happens if:
* The CIGAR is malformed, or doesn't match the sequence length, see [`cigar`](@ref).
* The auxiliary data is malformed or contains invalid values.
* The result is too long, e.g. a BAM record longer than 4 GiB.
* When converting to BAM, the reference or next reference of `src` is not present
  in the references of `dst` (`XAMErrors.unknown_reference`).

When converting to BAM, CIGARs with more than 65535 operations are stored in the
`CG` auxiliary tag, replacing any existing `CG` tag. When converting from BAM, a `CG` tag
holding the CIGAR is removed. Like htslib, sequence characters that can't be represented
in BAM, such as `.`, are converted to `N`.
"""
function try_copy! end

function try_copy!(dst::BAMRecord, src::SAMRecord)::Union{BAMRecord, XAMError}
    result = try_copy_kind!(dst, src)
    result isa BAMRecord && return result
    # `dst` may be partially overwritten, breaking its invariants
    reset!(dst)
    return result isa XAMErrorKind ? XAMError(result) : result
end

function try_copy_kind!(dst::BAMRecord, src::SAMRecord)::Union{BAMRecord, XAMError, XAMErrorKind}
    # Check what can fail before modifying `dst`
    ref_id = sam_reference_id(dst.refs, reference_name(src))
    ref_id isa XAMErrorKind && return ref_id
    next_ref_id = if is_rnext_equal(src)
        ref_id
    else
        sam_reference_id(dst.refs, next_reference_name(src))
    end
    next_ref_id isa XAMErrorKind && return next_ref_id
    cig = cigar(src)
    cig isa XAMError && return cig
    # `isvalid` also checks that the auxiliary data is well-formed
    src_aux = aux_data_unchecked(src)
    # Copying invalid auxiliary data throws
    isvalid(src_aux) || return XAMErrors.invalid_aux_data

    seq = seq_view(src)
    l_seq = is_missing(seq) ? 0 : length(seq)
    ref_len = cig === nothing ? 0 : CIGARStrings.ref_length(cig)
    # CIGARs with more than 65535 operations are stored in the CG tag. The CIGAR field
    # holds the placeholder `kSmN`, where k is the sequence length, and m the reference
    # length of the alignment, like htslib does. Each must fit in 28 bits.
    cigar_in_aux = cig !== nothing && length(cig) > typemax(UInt16)
    if cigar_in_aux && max(l_seq, ref_len) > 0x0fff_ffff
        return XAMErrors.record_too_long
    end
    n_cigar_op = if cig === nothing
        0
    elseif cigar_in_aux
        2
    else
        length(cig)
    end

    # Fixed-length fields. Positions are 0 if missing in SAM, -1 in BAM.
    dst.refID = ref_id
    dst.next_refID = next_ref_id
    dst.pos = (src.pos - 0x01) % Int32
    dst.next_pos = (src.pnext - 0x01) % Int32
    dst.mapq = src.mapq
    dst.flags = src.flags
    dst.tlen = src.tlen
    dst.n_cigar_op = n_cigar_op % UInt16
    dst.l_seq = l_seq % UInt32

    # Variable-length fields, except auxiliary data
    l_read_name = (src.qname_len % Int) + 1
    n_seq_bytes = cld(l_seq, 2)
    data = resize!(dst.data, l_read_name + 4 * n_cigar_op + n_seq_bytes + l_seq)
    mem = MemoryView(data)

    # Read name. A missing name is `*` in both SAM and BAM.
    copyto!(mem, @inbounds ImmutableMemoryView(src.data)[1:(src.qname_len % Int)])
    mem[l_read_name] = 0x00
    dst.l_read_name = l_read_name % UInt8
    offset = l_read_name

    # CIGAR
    if cigar_in_aux
        inbounds_store_le!(mem, offset + 1, ((l_seq % UInt32) << 4) | 0x04) # S
        inbounds_store_le!(mem, offset + 5, ((ref_len % UInt32) << 4) | 0x03) # N
    elseif cig !== nothing
        CIGARStrings.encode!(mem[(offset + 1):(offset + 4 * n_cigar_op)], BAMCIGAR, cig)
    end
    offset += 4 * n_cigar_op

    # Sequence and quality
    if l_seq > 0
        bamseq!(mem[(offset + 1):(offset + n_seq_bytes)], seq)
        offset += n_seq_bytes
        qual = quality_scores(src)
        if qual === nothing
            # Missing quality is encoded as all 0xff in BAM
            fill!(mem[(offset + 1):(offset + l_seq)], 0xff)
        else
            # Invariant: QUAL has the same length as SEQ
            copyto!(mem, offset + 1, qual, 1, l_seq)
        end
    end

    # Auxiliary data
    dst_aux = XAMAuxData.BAM.Auxiliary(data, length(data) + 1)
    copy!(dst_aux, src_aux)
    if cigar_in_aux
        cig = something(cig)
        delete!(dst_aux, "CG")
        # CG:B:I array of the CIGAR's BAM encoding
        CIGARStrings.encode!(append_cg_header!(data, length(cig)), BAMCIGAR, cig)
    end

    # Maintain the invariant that the block size fits in a UInt32
    length(data) > typemax(UInt32) - BAM_FIXED_BYTES && return XAMErrors.record_too_long
    return dst
end

function Base.copy!(dst::BAMRecord, src::SAMRecord)
    result = try_copy!(dst, src)
    result isa XAMError && throw(result)
    return result
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
        append!(data, codeunits(s))
        push!(data, UInt8('\t'))
        return length(data)
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

function try_copy!(dst::SAMRecord, src::BAMRecord)::Union{SAMRecord, XAMError}
    result = try_copy_kind!(dst, src)
    result isa SAMRecord && return result
    # `dst` may be partially overwritten, breaking its invariants
    reset!(dst)
    return result isa XAMErrorKind ? XAMError(result) : result
end

function Base.copy!(dst::SAMRecord, src::BAMRecord)
    result = try_copy!(dst, src)
    result isa XAMError && throw(result)
    return result
end

# Check if `cig`, the CIGAR of `record`, was loaded from the CG tag (see `load_aux_cigar`),
# as opposed to the CIGAR field.
function is_cigar_from_aux(record::BAMRecord, cig::BAMCIGAR)::Bool
    mem = MemoryView(cig)
    return GC.@preserve mem record begin
        pointer(mem) != pointer(record.data, cigar_offset(record) + 1)
    end
end

function try_copy_kind!(dst::SAMRecord, src::BAMRecord)::Union{SAMRecord, XAMError, XAMErrorKind}
    # Check what can fail before modifying `dst`
    cig = cigar(src)
    cig isa XAMError && return cig
    # `isvalid` also checks that the auxiliary data is well-formed
    src_aux = aux_data_unchecked(src)
    # Copying invalid auxiliary data throws
    isvalid(src_aux) || return XAMErrors.invalid_aux_data

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
    # The length includes the tab
    dst.rname_span = (rname_start, (rname_start + rname_len_p1 - 2) % UInt32)

    # POS. Invariant: src.pos is in BAM_POSITIONS, so this fits in SAM's POS
    dst.pos = src.pos + UInt32(1)
    append_number_field!(data, dst.pos, 10)

    # MAPQ
    append_number_field!(data, src.mapq, 10)
    dst.mapq = src.mapq

    # CIGAR
    dst.cigar_offset = length(data) % UInt32
    if cig === nothing
        push!(data, UInt8('*'))
    else
        CIGARStrings.encode_append!(data, CIGAR, cig)
    end
    push!(data, UInt8('\t'))

    # RNEXT
    if src.refID ≥ 0 && src.refID == src.next_refID
        # If RNAME and RNEXT is the same, we write it as `=`
        idx = (length(data) + 1) % UInt32
        dst.rnext_span = (idx, idx)
        push!(data, UInt8('='))
        push!(data, UInt8('\t'))
    else
        rnext_start = (length(data) + 1) % UInt32
        rnext_name_len_p1 = append_optional_string_field!(data, next_reference_name(src))
        dst.rnext_span = (rnext_start, (rnext_start + rnext_name_len_p1 - 2) % UInt32)
    end

    # PNEXT. Invariant: src.next_pos is in BAM_POSITIONS, so this fits in SAM's PNEXT
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
    qual = quality_string(src)
    if qual === nothing
        push!(data, UInt8('*'))
    else
        # `append!` would push each byte
        old_len = length(data)
        resize!(data, old_len + ncodeunits(qual))
        copyto!(data, old_len + 1, codeunits(qual), 1, ncodeunits(qual))
    end

    # Auxiliary fields
    dst.aux_offset = (length(data) + 1) % UInt32 # compensate for tab
    if !isempty(src_aux)
        push!(data, UInt8('\t'))
        dst_aux = XAMAuxData.SAM.Auxiliary(data, length(data) + 1)
        copy!(dst_aux, src_aux)
        # The CIGAR field holds the CIGAR, so it's not also stored in the CG tag
        if cig isa BAMCIGAR && is_cigar_from_aux(src, cig)
            delete!(dst_aux, "CG")
        end
        # If the CG tag was the only tag, remove the tab before it
        if isempty(dst_aux)
            pop!(data)
        end
    end

    # Invariant of `SAMRecord`
    length(data) ≤ MAX_SAM_LINE_LENGTH || return XAMErrors.line_too_long
    return dst
end
