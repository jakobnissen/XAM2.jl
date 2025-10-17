module XAM2

# TODO : Reset reader. Go back to after header, reset underlying.

using BufferIO: BufferIO,
    AbstractBufReader,
    BufReader,
    CursorReader,
    get_nonempty_buffer,
    consume,
    fill_buffer,
    get_buffer,
    line_views

using MemoryViews: ImmutableMemoryView, MemoryView, split_last, split_each, split_at
using StringViews: StringView
using BGZFLib: BGZFReader
using CIGARStrings: CIGARStrings, BAMCIGAR, CIGAR, CIGARError, unsafe_switch_memory

using XAMAuxData: XAMAuxData, AuxTag, try_auxtag

export BAMRecord, BAMReader

module XAMErrors

    @enum XAMError::UInt8 begin
        truncated_file
        invalid_block_length
        closed_stream
        header_truncated
        header_bad_magic_bytes
        header_top_bit_set
        invalid_cigar
    end
end # module XAMErrors

using .XAMErrors: XAMError

include("flags.jl")
include("stringvec.jl")
include("header.jl")
include("bamseq.jl")

struct RefMetaData end

struct Uninit end

# R can be reader or writer.
mutable struct BAMRecord{R, I}
    # This data contains, in order, with no unused bytes at the end:
    # 1. Read name including trailing null byte
    # 2. Cigar ops encoded in 32 bits each
    # 3. Seq encoded in 4 bits each
    # 4. Quality
    # 5. Auxiliary data
    const data::Vector{UInt8}
    # We need a ref here because the refids point back to the header, stored in the reader/writer
    const reader::R
    cigar::Union{Uninit, Nothing, BAMCIGAR, CIGARError}

    # fixed-length fields (see BAM specs for the details)
    # Length of alignment record, excluding this 32-bit field
    block_size::UInt32
    # Zero-based index of reference in , -1 if missing
    refid::Int32
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
    # Same as refid and pos
    next_refid::Int32
    next_pos::Int32
    # 0 if missing, else TLEN as defined in the SAM spec
    tlen::Int32

    global function new_uninitialized_record(::Type{I}, ::Type{R}, reader, data::Vector{UInt8}) where {I, R}
        return new{R, I}(data, reader, Uninit())
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
]

const BAM_FIXED_BYTES = length(BAM_RECORD_CONTENT)

@static if isdefined(Base, :memoryindex)
    memory_index(x::MemoryRef) = memoryindex(x)
else
    memory_index(x::MemoryRef) = Core.memoryrefoffset(x)
end

memory_index(v::Vector) = memory_index(Base.cconvert(Ptr, v))
memory_index(m::MemoryView) = first(only(parentindices(m)))

copy_cigar(x::Union{Nothing, Uninit, CIGARError}, _::Vector{UInt8}, _::Vector{UInt8}) = x

function copy_cigar(old_cigar::BAMCIGAR, old_v::Vector{UInt8}, cp::Vector{UInt8})
    # Find the view of `cp` that corresponds to the same data as the
    # view that backs `old_cigar`
    cigarmem = MemoryView(old_cigar)
    cigar_offset = memory_index(cigarmem)
    v_offset = memory_index(old_v)
    cigar_v_offset = cigar_offset - v_offset
    new_mem = ImmutableMemoryView(cp)[(cigar_v_offset + 1):(cigar_v_offset + length(cigarmem))]
    return unsafe_switch_memory(old_cigar, new_mem)
end

# Short form: Mapped, read name, strand
function Base.show(io::IO, record::BAMRecord)
    mapped = is_mapped(record) ? "Mapped" : "Unmapped"
    segment = is_forward(record) ? "Fw" : (is_reverse(record) ? "Rv" : "Md")
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

mutable struct BAMReader{T <: AbstractBufReader}
    const reader::BGZFReader{T}
    # Lens and names in their own fields, because these are accessed far more often
    # than any other metadata
    const ref_lens::Vector{Int32}
    const ref_names::StringVec
    const header::SAMHeader
    const offset_first_record::Int
    const copy::Bool
    error::Union{Nothing, XAMError}

    # These must be the last field so we can use a hack to initialize
    # the record and the reader at the same time.
    # That's also why it can't be const
    record::BAMRecord{BAMReader, T}

    global function new_bam_reader(reader::BGZFReader{T}, copy::Bool) where {T}
        buffer = UInt8[]
        # TODO: Handle potential errors - is maybe not SAMHEader, can be error
        # also, can throw
        (; header, stringvec, reflens, n_total_read) = parse_bam_header!(reader, buffer)

        result = new{T}(
            reader,
            reflens,
            stringvec,
            header,
            n_total_read,
            copy,
            nothing
        )
        result.record = new_bam_record(result, buffer)
        return result
    end
end

const FIXED_OFFSET = fieldoffset(BAMRecord{BAMReader, BGZFReader{CursorReader}}, 4)

function reset!(record::BAMRecord)
    resize!(record.data, 2)
    @inbounds record.data[1] = UInt8('*')
    @inbounds record.data[2] = 0x00
    record.cigar = Uninit()
    return reset_fixed_fields!(record)
end

@inline function reset_fixed_fields!(record::BAMRecord)
    GC.@preserve BAM_RECORD_CONTENT record begin
        dst = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            pointer(BAM_RECORD_CONTENT)::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return record
end

function Base.copy(record::BAMRecord{R, I}) where {R, I}
    cp = new_uninitialized_record(I, R, record.reader, copy(record.data))
    GC.@preserve cp record begin
        dst = Ptr{UInt8}(pointer_from_objref(cp)) + FIXED_OFFSET
        src = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            src::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return record
end

# TODO: Make this work!
function seek_first_record(reader::BAMReader)
    # Issue is that we need a virtual seek here, since we may not be in the beginning
    # of a block.
    seek(reader.reader, reader.offset_first_record)
    reset!(reader.record)
    return reader
end

function Base.iterate(reader::BAMReader, _::Nothing = nothing)
    result = try_read!(reader)
    result === nothing && return nothing
    isa(result, XAMError) && throw(result)
    record = reader.copy ? copy(result) : result
    return (record, nothing)
end

@inline function get_reader(x::BAMRecord{BAMReader, I}) where {I}
    return x.reader::BAMReader{I}
end

function parse_bam_header!(
        reader::BGZFReader,
        buffer::Vector{UInt8}
    )::Union{
        XAMError, @NamedTuple{
            header::SAMHeader,
            stringvec::StringVec,
            reflens::Vector{Int32},
            n_total_read::Int,
        }
    }

    memview = fill_bytes_to_data!(reader, buffer, 8)
    isnothing(memview) && return XAMErrors.header_truncated

    magic = unsafe_load_integer_le(UInt32, memview, 1)
    magic == 0x014d4142 || return XAMErrors.header_bad_magic_bytes

    l_text = unsafe_load_integer_le(UInt32, memview, 5)
    l_text > typemax(Int32) % UInt32 && return XAMErrors.header_top_bit_set
    n_total_read = 8

    memview = fill_bytes_to_data!(reader, buffer, l_text % Int + 4)
    isnothing(memview) && return XAMErrors.header_truncated

    # Copy the header because we need to store the data in the reader
    header_view = copy(@inbounds(memview[1:(end - 4)]))

    # TODO: This can throw. Is that alright? How do we handle it?
    header = parse_sam_header(header_view)
    stringvec = StringVec()
    reflens = Int32[]

    n_ref = unsafe_load_integer_le(UInt32, memview, l_text % Int + 1)
    n_ref > typemax(Int32) % UInt32 && return XAMErrors.header_top_bit_set
    n_total_read += (l_text % Int) + 4

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

        n_total_read += 4 + 4 + l_name
        push!(reflens, l_ref)
    end

    return (; header, stringvec, reflens, n_total_read)
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

function BAMReader(io::IO; copy::Bool = true, n_workers::Int = Threads.nthreads())
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
    record = new_uninitialized_record(T, BAMReader, reader, v)
    return reset_fixed_fields!(record)
end

flags(x::BAMRecord) = x.flags

is_mapped(x::BAMRecord) = !in(UNMAP, flags(x))

is_forward(x::BAMRecord) = in(READ1, flags(x))
is_reverse(x::BAMRecord) = in(READ2, flags(x))

function query_name(record::BAMRecord)::Union{AbstractString, Nothing}
    view = ImmutableMemoryView(record.data)
    view = view[1:Int(record.l_read_name - 0x01)]
    return if length(view) == 1 && @inbounds view[1] == UInt8('*')
        nothing
    else
        StringView(view)
    end
end

# TODO: Does this make sense? We already have ref information elsewhere.
# maybe delete it.
# TODO: Add more information?
# TODO: Add a show method
# This should be super cheap to construct - make it cheaper!
struct Reference
    name::StringView{ImmutableMemoryView{UInt8}}
    len::Int
end

function reference(record::BAMRecord)
    record.refid < 0 && return nothing
    id = record.refid + 1
    reader = get_reader(record)
    # TODO: A single boundscheck here?
    return Reference(
        reader.ref_names[id],
        reader.ref_lens[id]
    )
end

function reference_length(record::BAMRecord)::Union{Int, Nothing}
    record.refid < 0 && return nothing
    id = record.refid + 1
    reader = get_reader(record)
    return reader.ref_lens[id] % Int
end

struct HasEq end

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
    try_sequence(::Type{LongDNA{4}}, rec::BAMRecord)::Union{Nothing, HasEq, LongDNA{4}}

Get the sequence of `rec` as a `LongDNA{4}` biosequence.

This function returns `nothing` if the sequence is missing (i.e. if `sequence(rec)`
would return `nothing`)

This function returns `XAM2.HasEq()` if the sequence contains a `=`. In this case,
use the method `sequence(::LongDNA{4}, ::BAMRecord, ::BioSequence{<:DNAAlphabet})`.

!!! note
    This function requires `BioSequences` to be loaded
"""
function try_sequence end

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

function quality(record::BAMRecord)::Union{Nothing, ImmutableMemoryView{UInt8}}
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

function mapping_span(record::BAMRecord)
    left = @something left_position(record) return nothing
    cig = cigar(record)
    if cig === nothing || cig isa CIGARError
        return nothing
    end
    return left:(left + CIGARStrings.ref_length(cig) - 1)
end

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError}
    record.cigar === Uninit() || return record.cigar
    encoding = get_cigar_encoding(record)
    if encoding === nothing
        record.cigar = nothing
        return nothing
    end
    c = CIGARStrings.try_parse(BAMCIGAR, encoding)
    record.cigar = c
    return c
end

@inline function is_missing_encoding(record::BAMRecord, cigar_encoding::ImmutableMemoryView{UInt8})
    length(cigar_encoding) == 8 || return false
    u1 = unsafe_load_integer_le(UInt32, cigar_encoding, 1)
    u2 = unsafe_load_integer_le(UInt32, cigar_encoding, 5)
    return if ((u1 & 0x0f) != 0x04) | ((u2 & 0x0f) != 0x03)
        return false
    else
        ref_length = get_reader(record).ref_lens[record.refid + 1]
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

@inline function unsafe_load_integer_le(::Type{T}, mem::ImmutableMemoryView{UInt8}, i::Int) where {T <: Integer}
    x = zero(T)
    shift = 0
    mask = 8 * sizeof(T) - 1
    for offset in 0:(sizeof(T) - 1)
        byte = @inbounds mem[i + offset]
        x |= (byte % T) << (shift & mask)
        shift += 8
    end
    return x
end

function unsafe_load_fixed_fields(record::BAMRecord, buffer::ImmutableMemoryView{UInt8})
    GC.@preserve record buffer begin
        dst = Ptr{UInt8}(pointer_from_objref(record)) + FIXED_OFFSET
        @ccall memcpy(
            dst::Ptr{Cchar},
            pointer(buffer)::Ptr{Cchar},
            (BAM_FIXED_BYTES % UInt)::Csize_t,
        )::Cvoid
    end
    return nothing
end

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

function reg2bin(span::UnitRange{Int32})::UInt16
    # TODO: remove these checks, also the one at the bottom
    (start, stop) = (first(span), last(span))
    @assert start in Int32(0):Int32(2^29 - 1)
    @assert stop in start:(Int32(2^29 - 1))

    # If start is zero, the read is unmapped, and so the range
    # must be 0:0 exactly
    @assert !iszero(start) || iszero(stop)

    # Convert start to an offset to make it zero-indexed.
    # Note that this will be -1 when span is 0:0.
    # Then this will use arithmetic shit right.
    start -= one(start)
    result = if start >> 14 == stop >> 14
        div((1 << 15) - 1, 7) + (start >> 14)
    elseif start >> 17 == stop >> 17
        div((1 << 12) - 1, 7) + (start >> 17)
    elseif start >> 20 == stop >> 20
        div((1 << 9) - 1, 7) + (start >> 20)
    elseif start >> 23 == stop >> 23
        div((1 << 6) - 1, 7) + (start >> 23)
    elseif start >> 26 == stop >> 26
        div((1 << 3) - 1, 7) + (start >> 26)
    else
        0
    end

    @assert result in 0:Int(typemax(UInt16))

    return result % UInt16
end

function reg2bins!(bins::Vector{UInt16}, span::UnitRange{Int32})
    (start, stop) = (first(span) - Int32(1), last(span) - Int32(1))
    # @assert start in Int32(1):Int32(2^29 - 1)
    # @assert stop in start:(Int32(2^29 - 1))

    # If start is zero, the read is unmapped, and so the range
    # must be 0:0 exactly
    # @assert !iszero(start) || iszero(stop)

    # TODO: Check the 0,0 case

    ranges = map((1, 9, 73, 585, 4681)) do idx
        idx = idx % UInt16
        (idx + (start >> 26) % UInt16):(idx + (stop >> 26) % UInt16)
    end

    resize!(bins, sum(length, ranges; init = 0) + 1)

    @inbounds bins[1] = 0

    i = 1
    for range in ranges, k in range
        i += 1
        @inbounds bins[i] = k % UInt16
    end

    return bins
end

function tryparse_base10(::Type{T}, mem::ImmutableMemoryView{UInt8}) where {T <: Union{UInt8, UInt16, UInt32, UInt64}}
    isempty(mem) && return nothing
    n = UInt64(0)
    has_overflowed = false
    bad_byte = false
    for i in mem
        bad_byte |= i ∉ 0x30:0x39
        (n, f) = Base.Checked.mul_with_overflow(n, UInt64(0x0a))
        has_overflowed |= f
        (n, f) = Base.Checked.add_with_overflow(n, (i - 0x30) % UInt64)
        has_overflowed |= f
    end
    ((n > typemax(T)) | has_overflowed | bad_byte) && return nothing
    return n % T
end


end # module XAM2
