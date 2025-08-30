module XAM2

using BufIO: AbstractBufReader, BufReader, get_nonempty_buffer
using MemoryViews: ImmutableMemoryView, MemoryView, split_last
using StringViews: StringView
using BGZFLib: BGZFReader
using BufIO: CursorReader, consume, fill_buffer, get_buffer
using CIGARStrings: CIGARStrings, BAMCIGAR, CIGAR, CIGARError
using XAMAuxData: XAMAuxData

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

const BAM_FIXED_BYTES = 36

struct RefMetaData end

# R can be reader or writer.
mutable struct BAMRecord{R, I}
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
    # This data contains, in order, with no unused bytes at the end:
    # 1. Read name including trailing null byte
    # 2. Cigar ops encoded in 32 bits each
    # 3. Seq encoded in 4 bits each
    # 4. Quality
    # 5. Auxiliary data
    const data::Vector{UInt8}
    # We need a ref here because the refids point back to the header, stored in the reader/writer
    const reader::R
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
    refname = reference_name(record)
    if refname === nothing
        println(buf, "< missing >")
    else
        println(buf, rtruncate(repr(refname), width - 12))
    end

    # Mapping span
    print(buf, "  Ref pos: ")
    left = left_position(record)
    right = right_position(record)
    if right !== nothing
        println(buf, left, ':', right)
    elseif left !== nothing
        println(buf, left)
    else
        println(buf, "< missing >")
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
    pop!(v)
    return write(io, v)
end

mutable struct BAMReader{T <: IO}
    reader::BGZFReader{T}
    # Lens and names in their own fields, because these are accessed far more often
    # than any other metadata
    ref_lens::Vector{Int32}
    ref_names::StringVec
    ref_metadata::RefMetaData
    copy::Bool
    error::Union{Nothing, XAMError}

    # This must be the last field so we can use a hack to initialize
    # the record and the reader at the same time
    record::BAMRecord{BAMReader, T}

    global function new_bam_reader(reader::BGZFReader{T}, copy::Bool) where {T}
        result = new{T}(
            reader,
            Int32[],
            StringVec(),
            RefMetaData(),
            copy,
            nothing
        )
        record = BAMRecord(result)
        result.record = record
        return result
    end
end

@inline function get_reader(x::BAMRecord{BAMReader, I}) where I
    return x.reader::BAMReader{I}
end

# TODO: Handle errors better
function init_bam_reader!(reader::BAMReader)
    memview = fill_bytes_to_data!(reader, 8)
    if isnothing(memview)
        reader.error = XAMErrors.header_truncated
        return reader
    end

    magic = unsafe_load_integer_le(UInt32, memview, 1)
    if magic != 0x014d4142
        reader.error = XAMError.header_bad_magic_bytes
        return reader
    end

    l_text = unsafe_load_integer_le(UInt32, memview, 5)
    if l_text > typemax(Int32) % UInt32
        reader.error = XAMError.header_top_bit_set
        return reader
    end

    memview = fill_bytes_to_data!(reader, l_text % Int + 4)
    if isnothing(memview)
        reader.error = XAMErrors.header_truncated
        return reader
    end

    # TODO: Parse header - we skip for now

    n_ref = unsafe_load_integer_le(UInt32, memview, length(memview) - 3)
    if n_ref > typemax(Int32) % UInt32
        reader.error = XAMError.header_top_bit_set
        return reader
    end

    for _ in 1:(n_ref % Int)
        memview = fill_bytes_to_data!(reader, 4)
        if isnothing(memview)
            reader.error = XAMErrors.header_truncated
            return reader
        end

        l_name = unsafe_load_integer_le(UInt32, memview, 1)

        memview = fill_bytes_to_data!(reader, l_name + 4)
        if isnothing(memview)
            reader.error = XAMErrors.header_truncated
            return reader
        end

        # Skip null byte
        push_codeunits!(reader.ref_names, memview[1:(l_name - 1)])

        l_ref = unsafe_load_integer_le(Int32, memview, l_name + 1)

        if l_ref < 1
            reader.error = XAMErrors.header_empty_ref
            return reader
        end

        push!(reader.ref_lens, l_ref)
    end
    return reader
end

function fill_bytes_to_data!(reader::BAMReader, nbytes::Int)::Union{ImmutableMemoryView{UInt8}, Nothing}
    data = reader.record.data
    length(data) < nbytes && resize!(data, nbytes)
    remaining = MemoryView(data)[1:nbytes]
    while !isempty(remaining)
        buffer = @something get_nonempty_buffer(reader.reader) return nothing
        mn = min(length(buffer), length(remaining))
        @assert mn > 0
        @assert mn <= length(buffer)
        slice = @inbounds buffer[1:mn]
        copyto!(remaining, slice)
        consume(reader.reader, mn)
        remaining = @inbounds remaining[(mn + 1):end]
    end
    return @inbounds ImmutableMemoryView(data)[1:nbytes]
end

Base.isopen(x::BAMReader) = isopen(x.reader)

function BAMReader(reader::BGZFReader{T}; copy::Bool = true) where {T}
    return init_bam_reader!(new_bam_reader(reader, copy))
end

function BAMReader(io::IO; copy::Bool = true, threads::Int = Threads.nthreads())
    reader = BGZFReader(io; threads)
    return init_bam_reader!(new_bam_reader(reader, copy))
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

function BAMRecord(reader::BAMReader{T}) where T
    return BAMRecord{BAMReader, T}(
        BAM_FIXED_BYTES - sizeof(Int32) + 2, # size of empty BAM record
        -1, # meaning unmapped
        -1, # unmapped
        2, # read name '*' + null byte
        0xff, # mapping quality (unavailable)
        4680, # TODO: Understand why, but BAM specs uses this value
        0, # cigar ops
        new_sam_flags(0x0004), # Flag: Unmapped
        0, # length of seq
        -1,
        -1,
        0,
        UInt8[UInt8('*'), 0x00],
        reader,
    )
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

function reference_name(record::BAMRecord)::Union{AbstractString, Nothing}
    record.refid < 0 && return nothing
    return get_reader(record).ref_names[record.refid + 1]
end

struct HasEq end

function has_eq_seq_symbol(record::BAMRecord)
    iszero(record.l_seq) && return false
    n_seq_whole_bytes = record.l_seq >> 1
    offset = record.l_read_name + 4 * record.n_cigar_op
    has_half_byte = isodd(record.l_seq)
    span = (offset + 1):(offset + n_seq_whole_bytes + has_half_byte)
    mem = ImmutableMemoryView(record.data)[span]
    if has_half_byte
        (last_byte, mem) = @inbounds split_last(mem)
        iszero(last_byte & 0xf0) && return true
    end
    i = 1
    chunk_size = 128
    @inbounds while i <= length(mem) - chunk_size + 1
        any = false
        for j in i:(i + chunk_size - 1)
            any |= iszero(mem[j] & 0x0f) | iszero(mem[j] & 0xf0)
        end
        any && return true
        i += chunk_size
    end
    @inbounds while i <= length(mem)
        (iszero(mem[i] & 0x0f) | iszero(mem[i] & 0xf0)) && return true
        i += 1
    end
    return false
end

function try_sequence end

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

function right_position(record::BAMRecord)::Union{Int, Nothing}
    left = @something left_position(record) return nothing
    cig = cigar(record)
    if cig === nothing || cig isa CIGARError
        return nothing
    end
    return left + CIGARStrings.ref_length(cig) - 1
end

function cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError}
    encoding = @something get_cigar_encoding(record) return nothing
    return CIGARStrings.try_parse(BAMCIGAR, encoding)
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
    view = @inbounds view[(x.l_read_name + 1):lastindex]
    return is_missing_encoding(x, view) ? load_aux_cigar(x) : view
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
    GC.@preserve record buffer unsafe_copyto!(
        Ptr{UInt8}(pointer_from_objref(record)), pointer(buffer), BAM_FIXED_BYTES
    )
    nothing
end

function try_read!(x::BAMReader)::Union{BAMRecord, Nothing, XAMError}
    isnothing(x.error) || return x.error

    # Check if EOF and return that error if so
    buffer = get_nonempty_buffer(x.reader)
    isnothing(buffer) && return nothing

    record = x.record
    vector = record.data

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


end # module XAM2
