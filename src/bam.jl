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
as an error value, e.g. [`cigar`](@ref) returns a `CIGARError` for a malformed CIGAR.

Auxiliary data is returned as an `XAMAuxData.BAM.Auxiliary`, which may throw
when accessed if it is malformed. Use `XAMAuxData.is_well_formed` to check this first.

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

    # This field is cached. Uninit means: Not computed. Nothing means no cigar.
    # Otherwise, a cigar or CIGARError for malformed CIGAR, or XAMError for CIGAR
    # whose length does not match the query length.
    cigar::Union{Uninit, Nothing, BAMCIGAR, CIGARError, XAMError}

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

    global function new_uninitialized_record(refs::References, data::Vector{UInt8})
        return new(data, refs, Uninit())
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
If the binary header is malformed, `io` is closed and an `XAMError` is thrown.
See [`XAM2.try_init`](@ref) for a version that returns the error instead.
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
    buffer = UInt8[]
    parsed = close_on_throw(() -> parse_bam_header!(reader, buffer), reader)
    if parsed isa XAMError
        close(reader)
        return parsed
    end
    (; header, refs) = parsed
    if strict_header && header isa XAMError
        close(reader)
        return header
    end
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
    # The BGZFReader constructor reads from `io`, which may throw
    bgzf = close_on_throw(() -> BGZFReader(io; n_workers), io)
    return try_init(BAMReader, bgzf; copy, strict_header)
end

get_header(reader::BAMReader) = reader.header
references(reader::BAMReader) = reader.refs

const FIXED_BAM_OFFSET = fieldoffset(BAMRecord, 4)

@assert sizeof(BAMRecord) == FIXED_BAM_OFFSET + BAM_FIXED_BYTES

"""
    reset!(::BAMRecord)

Resets the `BAMRecord` to its default state, deleting all information from it.
After this function, the record is equal to an empty record created with
`BAMRecord(record.refs)`.
"""
function reset!(record::BAMRecord)
    resize!(record.data, 2)
    @inbounds record.data[1] = UInt8('*')
    @inbounds record.data[2] = 0x00
    record.cigar = Uninit()
    return reset_fixed_fields!(record)
end

Base.IteratorSize(::Type{<:BAMReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:BAMReader}) = BAMRecord

function Base.iterate(reader::BAMReader, _::Nothing = nothing)
    result = try_read!(reader)
    result === nothing && return nothing
    isa(result, XAMError) && throw(result)
    record = reader.copy ? copy(result) : result
    return (record, nothing)
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
    vector_copy = copy(record.data)
    cp = new_uninitialized_record(record.refs, vector_copy)
    GC.@preserve cp record begin
        dst = Ptr{UInt8}(pointer_from_objref(cp)) + FIXED_BAM_OFFSET
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_BAM_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            src::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    cp.cigar = copy_cigar!(record.cigar, record.data, vector_copy)
    return cp
end

function seek_first_record(reader::BAMReader)
    return seek(reader, reader.offset_first_record)
end

function Base.seek(reader::BAMReader, offset::VirtualOffset)
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

Base.isopen(x::BAMReader) = isopen(x.reader)

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


function BAMReader(f::Function, args...; kwargs...)
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

function Base.close(reader::BAMReader)
    reader.error = XAMError(XAMErrors.closed_stream)
    return close(reader.reader)
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
    sequence(rec::BAMRecord)::Union{Nothing, AbstractString}

Get the query sequence of `rec` as an `AbstractString`, or
`nothing` if the sequence is not available.

The returned `AbstractString`'s type is an implementation detail, but it is guaranteed to:
* Be constructed cheaply in constant time
* Backed by a view into the existing data of the record (so not separately allocated)
* Contain at least 1 codeunit, and having the codeunit type `UInt8`.
* Containing only the `Char`s `"=ACMGRSVTWYHKDBN"` encoded in ASCII

!!! warning
    The returned sequence aliases the record. Mutating the record will
    invalidate the sequence
"""
function sequence(record::BAMRecord)::Union{BAMSeq, Nothing}
    iszero(record.l_seq) && return nothing
    span = (seq_offset(record) + 1):qual_offset(record)
    return BAMSeq(@inbounds(ImmutableMemoryView(record.data)[span]), record.l_seq % Int)
end

"""
    quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}

Return the PHRED quality scores of `record`, or `nothing` if the quality is
missing (stored as all `0xff`, or the sequence is missing).

The scores have no offset, so e.g. a PHRED score of 31 is `UInt8(31)`, as opposed to
`UInt8(31 + 33)` in SAM. They are stored like this in BAM, so this is a view of the
record's data.

See also: [`quality_string`](@ref)

!!! warning
    The returned view aliases the record. Mutating the record will
    invalidate the view.
"""
function quality_scores(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
    iszero(record.l_seq) && return nothing
    span = (qual_offset(record) + 1):aux_offset(record)
    qual = @inbounds ImmutableMemoryView(record.data)[span]

    # Missing quality is encoded as all-0xFF
    return if @inbounds first(qual) == 0xff
        nothing
    else
        qual
    end
end

"""
    quality_string(record::BAMRecord)::Union{Nothing, StringView{QualityView}}

Return the quality of `record` as a PHRED+33-encoded string, as in the QUAL field of
SAM, or `nothing` if the quality is missing (see [`quality_scores`](@ref)).

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
    aux_data(record::BAMRecord)::XAMAuxData.BAM.MutableAuxiliary

Return the auxiliary data as a mutable `XAMAuxData.BAM.Auxiliary`.
If the record does not contain any auxiliary fields, the returned object
is empty.

!!! warning
    The returned `Auxiliary` object aliases the record.
    Mutating the record will invalidate the auxiiary data.

See the documentation of `XAMAuxData.BAM.Auxiliary` for more details.
"""
function aux_data(record::BAMRecord)
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

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError, XAMError}
    record.cigar === Uninit() || return record.cigar
    encoding = get_cigar_encoding(record)
    if encoding === nothing
        record.cigar = nothing
        return nothing
    end
    cig = CIGARStrings.try_parse(BAMCIGAR, encoding)
    if !isa(cig, CIGARError)
        cig = check_query_length(cig, record.l_seq % Int)
    end
    return record.cigar = cig
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
    # Accessing malformed aux data may throw, so we check first.
    is_well_formed(aux) || return placeholder
    v = get(aux, "CG", nothing)
    return v isa AbstractVector{UInt32} ? parent(v) : placeholder
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
    record.cigar = Uninit()

    # Read the block size, and the fixed-length fields into the record
    memview = fill_bytes_to_data!(x, 4 + BAM_FIXED_BYTES)
    memview isa ImmutableMemoryView || return fail_read!(x, fill_error(memview, offset))
    stored_block_size = unsafe_load_integer_le(UInt32, memview, 1)
    unsafe_load_fixed_fields(record, @inbounds memview[5:end])
    # Clear undefined flag bits, which `SAMFlags` must not have
    record.flags = SAMFlags(Integer(record.flags))
    err = check_fixed_fields(record, stored_block_size)
    err === nothing || return fail_read!(x, XAMError(err; offset))

    # Read all the variable length fields, overwriting the fixed field bytes,
    # which we have already moved into the integer fields of the struct
    data_len = (stored_block_size - UInt32(BAM_FIXED_BYTES)) % Int
    memview = fill_bytes_to_data!(x, data_len)
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
    write(writer, htol(0x014d4142))
    write(writer, htol(length(header_text) % UInt32))
    write(writer, header_text)
    write(writer, htol(length(refs) % UInt32))

    # Write reference sequences
    for ref in refs
        cu = codeunits(ref.name)
        write(writer, htol((length(cu) + 1) % UInt32))
        write(writer, cu, 0x00)
        write(writer, htol(ref.length % UInt32))
    end

    maps = IdDict{References, Union{Nothing, Vector{Int32}}}()
    return BAMWriter{T}(writer, refs, nothing, nothing, maps)
end

function check_writer_header(header::SAMHeader, refs::References)::Nothing
    if length(MemoryView(header)) > typemax(Int32)
        throw(ArgumentError("Header must be 2147483647 bytes or shorter"))
    end
    check_header_references(header, refs)
    return nothing
end

function BAMWriter(f::Function, args...; kwargs...)
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

Base.close(writer::BAMWriter) = close(writer.writer)
Base.flush(writer::BAMWriter) = flush(writer.writer)

"""
    XAM2.try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMError}

Write `record` to `writer`, and return the number of bytes written.
If the record can't be represented using the references of `writer`, write nothing,
and return an `XAMError`. This happens if the reference or next reference of
`record` is not present in `writer`, or if it has a different length in `writer`.
Errors from the underlying IO are thrown.

See also: [`BAMWriter`](@ref)
"""
function try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMError}
    result = try_write_kind(writer, record)
    result isa XAMErrorKind && return XAMError(result)
    return result
end

function try_write_kind(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMErrorKind}
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
    return write_record(writer, record, ref_id, next_ref_id)
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

# Write the record, but with the given reference IDs instead of the record's own.
function write_record(writer::BAMWriter, record::BAMRecord, ref_id::Int32, next_ref_id::Int32)::Int
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
    # Byte offsets of refID and next_refID, after the block size
    inbounds_store_le!(buffer, 5, ref_id)
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
    next_ref_id = sam_reference_id(dst.refs, next_reference_name(src))
    next_ref_id isa XAMErrorKind && return next_ref_id
    cig = cigar(src)
    cig isa CIGARError && return XAMErrors.invalid_cigar
    cig isa XAMError && return cig
    src_aux = aux_data(src)
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
    # Like htslib, records with no CIGAR, or one spanning no reference positions,
    # are binned as spanning one position. Unmapped records with no position get bin 4680.
    dst.bin = if dst.pos < 0
        UInt16(4680)
    else
        # As `UInt64`, since the end may not fit in 32 bits
        beg = dst.pos % UInt64
        reg2bin(beg, beg + (max(ref_len, 1) % UInt64))
    end
    # The CIGAR is not cached, because it would alias `data`, which is resized below
    dst.cigar = Uninit()

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
        delete!(dst_aux, "CG")
        # CG:B:I array of the CIGAR's BAM encoding
        old_len = length(data)
        resize!(data, old_len + 8 + 4 * length(cig))
        mem = MemoryView(data)[(old_len + 1):end]
        copyto!(mem, b"CGBI")
        inbounds_store_le!(mem, 5, length(cig) % UInt32)
        CIGARStrings.encode!(mem[9:end], BAMCIGAR, cig)
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
    cig isa CIGARError && return XAMErrors.invalid_cigar
    cig isa XAMError && return cig
    src_aux = aux_data(src)
    # Copying invalid auxiliary data throws
    isvalid(src_aux) || return XAMErrors.invalid_aux_data

    data = dst.data
    empty!(data)
    # The CIGAR is not cached, because it would alias `data`, which is resized below
    dst.cigar = Uninit()

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
