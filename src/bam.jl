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
    # 5. `pos >= -1` and `next_pos >= -1`
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
    # Otherwise, a cigar or invalid data.
    cigar::Union{Uninit, Nothing, BAMCIGAR, CIGARError, XAMErrorKind}

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
    elseif cig isa Union{CIGARError, XAMErrorKind}
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
    const refs::References
    const header::SAMHeader
    const offset_first_record::VirtualOffset
    const copy::Bool
    const record::BAMRecord
    error::Union{Nothing, XAMErrorKind}
end

function new_bam_reader(reader::BGZFReader{T}, copy::Bool) where {T}
    buffer = UInt8[]
    # TODO: Add a non-throwing constructor (see REVIEW.md, section 4)
    parsed = parse_bam_header!(reader, buffer)
    parsed isa SAMHeaderError && throw(parsed)
    parsed isa XAMErrorKind && throw_xam_error(parsed)
    (; header, refs) = parsed
    vo_at_first_record = virtual_position(reader)
    record = new_bam_record(refs, buffer)
    return BAMReader{T}(reader, refs, header, vo_at_first_record, copy, record, nothing)
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
    isa(result, XAMErrorKind) && throw_xam_error(result)
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
    seek(reader.reader, reader.offset_first_record)
    reset!(reader.record)
    # Errors are sticky until the reader is seeked, except if it's closed
    reader.error == XAMErrors.closed_stream || (reader.error = nothing)
    return reader
end

function parse_bam_header!(
        reader::BGZFReader,
        buffer::Vector{UInt8}
    )::Union{
        SAMHeaderError, XAMErrorKind, @NamedTuple{
            header::SAMHeader,
            refs::References,
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

    header = try_parse_header(header_view)
    header isa SAMHeaderError && return header
    stringvec = StringVec()
    reflens = Int32[]

    n_ref = unsafe_load_integer_le(UInt32, memview, l_text % Int + 1)
    n_ref > typemax(Int32) % UInt32 && return XAMErrors.header_top_bit_set

    for _ in 1:(n_ref % Int)
        memview = fill_bytes_to_data!(reader, buffer, 4)
        isnothing(memview) && return XAMErrors.header_truncated

        # Name length includes the trailing null byte
        l_name = unsafe_load_integer_le(UInt32, memview, 1)
        in(l_name, 2:(typemax(Int32) % UInt32)) || return XAMErrors.invalid_reference
        l_name = l_name % Int

        memview = fill_bytes_to_data!(reader, buffer, l_name + 4)
        isnothing(memview) && return XAMErrors.header_truncated

        @inbounds(memview[l_name]) == 0x00 || return XAMErrors.invalid_reference
        push_codeunits!(stringvec, @inbounds memview[1:(l_name - 1)])

        l_ref = unsafe_load_integer_le(Int32, memview, l_name + 1)
        push!(reflens, l_ref)
    end

    refs = try_build_references(stringvec, reflens)
    refs isa XAMErrorKind && return refs
    return (; header, refs)
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
    reader.error = XAMErrors.closed_stream
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

function first_eq_seq_pos(record::BAMRecord)
    iszero(record.l_seq) && return nothing
    n_seq_whole_bytes = (record.l_seq >> 1) % Int
    offset = seq_offset(record)
    has_half_byte = isodd(record.l_seq)
    span = (offset + 1):(offset + n_seq_whole_bytes + has_half_byte)
    mem = @inbounds ImmutableMemoryView(record.data)[span]
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

Return the quality as a PHRED-encoded memory view, or `nothing` if the quality is
missing (set to all `0xFF`, or the sequence is missing).

Note that the returned quality does not have a PHRED offset, so e.g. a PHRED
score of 31 is encoded as `UInt8(31)` as opposed to e.g. `UInt8(31 + 33)`.

!!! warning
    The returned array aliases the record. Mutating the record will
    invalidate the array
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

function left_position(record::BAMRecord)::Union{Int, Nothing}
    result = record.pos + 1
    return iszero(result) ? nothing : result
end

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError, XAMErrorKind}
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
# reference length.
@inline function is_placeholder_cigar(record::BAMRecord, cigar_encoding::ImmutableMemoryView{UInt8})
    length(cigar_encoding) == 8 || return false
    record.refID < 0 && return false
    u1 = unsafe_load_integer_le(UInt32, cigar_encoding, 1)
    u2 = unsafe_load_integer_le(UInt32, cigar_encoding, 5)
    return if ((u1 & 0x0f) != 0x04) | ((u2 & 0x0f) != 0x03)
        false
    else
        ref_length = @inbounds record.refs.lengths[record.refID + 1]
        ((u1 >> 4) == record.l_seq) & ((u2 >> 4) == ref_length)
    end
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
    XAM2.try_read!(reader::BAMReader)::Union{BAMRecord, Nothing, XAMErrorKind}

Attempt to read the next record of `reader`, similar to iterating the reader.
If the reader is EOF, return `nothing`.
If the next record is malformed, return a `XAMErrorKind`.
Once this function returns a `XAMErrorKind`, subsequent calls to this function
will return the same error until the reader is closed or seeked.

The structure of the record is validated, e.g. that its fields fit in the record,
and that its reference IDs and positions are in bounds. The content of variable-length
fields, e.g. the CIGAR, is validated when it is accessed.

This function does not throw if the reader's data is malformed.
"""
function try_read!(x::BAMReader)::Union{BAMRecord, Nothing, XAMErrorKind}
    isnothing(x.error) || return x.error

    # Check if EOF and return that error if so
    buffer = get_nonempty_buffer(x.reader)
    isnothing(buffer) && return nothing

    record = x.record
    record.cigar = Uninit()

    # Read the block size, and the fixed-length fields into the record
    memview = fill_bytes_to_data!(x, 4 + BAM_FIXED_BYTES)
    isnothing(memview) && return fail_read!(x, XAMErrors.truncated_file)
    stored_block_size = unsafe_load_integer_le(UInt32, memview, 1)
    unsafe_load_fixed_fields(record, @inbounds memview[5:end])
    err = check_fixed_fields(record, stored_block_size)
    err === nothing || return fail_read!(x, err)

    # Read all the variable length fields, overwriting the fixed field bytes,
    # which we have already moved into the integer fields of the struct
    data_len = (stored_block_size - UInt32(BAM_FIXED_BYTES)) % Int
    memview = fill_bytes_to_data!(x, data_len)
    isnothing(memview) && return fail_read!(x, XAMErrors.truncated_file)
    resize!(record.data, data_len)
    @inbounds(record.data[record.l_read_name]) == 0x00 || return fail_read!(x, XAMErrors.invalid_read_name)
    return record
end

# Set the reader's sticky error, and reset its record, since a failed read
# may leave the record in a state that violates its invariants.
@noinline function fail_read!(x::BAMReader, err::XAMErrorKind)::XAMErrorKind
    reset!(x.record)
    x.error = err
    return err
end

# Check the structural invariants of `record` that only depend on its fixed fields,
# see the comment in the definition of `BAMRecord`.
function check_fixed_fields(record::BAMRecord, stored_block_size::UInt32)::Union{Nothing, XAMErrorKind}
    # The block size includes the fixed fields after the block size field.
    stored_block_size < BAM_FIXED_BYTES && return XAMErrors.invalid_block_length
    data_len = (stored_block_size - UInt32(BAM_FIXED_BYTES)) % Int
    iszero(record.l_read_name) && return XAMErrors.invalid_read_name
    n_refs = length(record.refs)
    is_valid_reference_id(record.refID, n_refs) || return XAMErrors.invalid_reference_id
    is_valid_reference_id(record.next_refID, n_refs) || return XAMErrors.invalid_reference_id
    (record.pos < -1 || record.next_pos < -1) && return XAMErrors.invalid_position
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
    BAMWriter(open("out.bam", "w"), get_header(reader), references(reader)) do writer
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
    XAM2.try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMErrorKind}

Write `record` to `writer`, and return the number of bytes written.
If the record can't be represented using the references of `writer`, write nothing,
and return an `XAMErrorKind`. This happens if the reference or next reference of
`record` is not present in `writer`, or if it has a different length in `writer`.
Errors from the underlying IO are thrown.

See also: [`BAMWriter`](@ref)
"""
function try_write(writer::BAMWriter, record::BAMRecord)::Union{Int, XAMErrorKind}
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
    result isa XAMErrorKind && throw_xam_error(result)
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

function get_index_of_refname(record::BAMRecord, name::Union{AbstractString, Nothing})
    name === nothing && return Int32(-1)
    id = reference_id(record.refs, name)
    # TODO: Return an error value instead
    id === nothing && throw_xam_error(XAMErrors.unknown_reference)
    return id
end

function Base.copy!(dst::BAMRecord, src::SAMRecord)
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
        if !isa(ms, UnitRange{Int})
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
    elseif cig isa CIGARError
        throw_xam_error(XAMErrors.invalid_cigar)
    elseif cig isa XAMErrorKind
        throw_xam_error(cig)
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
            reflen = dst.refs.lengths[dst.refID + 1]
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

    # Maintain the invariant that the block size fits in a UInt32
    if length(dst.data) > typemax(UInt32) - BAM_FIXED_BYTES
        reset!(dst)
        throw_xam_error(XAMErrors.record_too_long)
    end

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
        throw_xam_error(XAMErrors.invalid_cigar)
    elseif cig isa XAMErrorKind
        throw_xam_error(cig)
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
