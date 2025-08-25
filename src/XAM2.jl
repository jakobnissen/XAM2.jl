module XAM2

# TODO:
#=
* Use BGZFLib, and make reader take one of those
* Make a constructor for Reader, making the record/reader self-reference work
* Try reading a bunch of records

* Make a new type for the sequence, make it <: AbstractString
* In extension, make a constructor for 4-bit nucleotide for sequence

=#

using BufIO: AbstractBufReader, BufReader, get_nonempty_buffer
using BGZFStreams: BGZFStreams
using MemoryViews: ImmutableMemoryView, MemoryView
using StringViews: StringView

module XAMErrors

@enum XAMError::UInt8 EndOfFile TruncatedFile InvalidBlockLength

end # module XAMErrors

using .XAMErrors: XAMError

include("flags.jl")
include("stringvec.jl")

const BAM_FIXED_BYTES = 36

struct RefMetaData end # TODO

# R can be reader or writer.
# TODO: What do we do about sequence and CIGAR?
# Maybe keep a lazily initialized field for CIGAR,
# and a new BioSequence type with a special constructor for
# LongSequence (see XAM.jl)
mutable struct BAMRecord{R}
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

mutable struct BAMReader{R <: AbstractBufReader}
    reader::R
    # Lens and names in their own fields, because these are accessed far more often
    # than any other metadata
    ref_lens::Vector{Int32}
    ref_names::StringVec
    ref_metadata::RefMetaData
    copy::Bool

    # This must be the last field so we can use a hack to initialize
    # the record and the reader at the same time
    record::BAMRecord
end

function BAMReader(io::IO; copy::Bool = true)
    # TODO: use BGZFLib, also set threads in BAMReader
    bgzf = BGZFStream(io)
    reader = BufReader(bgzf)
    return BAMReader{typeof(reader)}(reader, copy)
end

function BAMRecord(reader::BAMReader)
    return BAMRecord{BAMReader}(
        BAM_FIXED_BYTES - sizeof(Int32) + 2, # size of empty BAM record
        -1, # meaning unmapped
        -1, # unmapped
        2, # read name 'A' + null byte
        0xff, # mapping quality (unavailable)
        4680, # TODO: Understand why, but BAM specs uses this value
        0, # cigar ops
        4, # Flag: Unmapped
        0, # length of seq
        -1,
        -1,
        0,
        UInt8[0x61, 0x00],
        reader,
    )
end

flags(x::BAMRecord) = x.flags

is_mapped(x::BAMRecord) = !in(UNMAP, flags(x))

# If more than typemax(UInt16) cigar ops, field is missing and the data is encoded
# in the aux
function is_missing_cigar_field(x::BAMRecord)
    x.n_cigar_op == 2 || return false
    is_mapped(x) || return false
    ref_length = x.reader.ref_lens[x.refid + 1]
    return (x.data[x.l_read_name + 1] == (x.l_seq % UInt32 << 4 | 4)) &
        (x.data[x.l_read_name + 2] == (ref_length % UInt32 << 4 | 3))
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
    x
end

function unsafe_load_fixed_fields(record::BAMRecord, buffer::ImmutableMemoryView)
    record.block_size = unsafe_load_integer_le(UInt32, buffer, 1)
    record.refid = unsafe_load_integer_le(Int32, buffer, 5)
    record.pos = unsafe_load_integer_le(Int32, buffer, 9)
    record.l_read_name = @inbounds buffer[13]
    record.mapq = @inbounds buffer[14]
    record.bin = unsafe_load_integer_le(UInt16, buffer, 15)
    record.n_cigar_op = unsafe_load_integer_le(UInt16, buffer, 17)
    record.flags = SAMFlags(unsafe_load_integer_le(UInt16, buffer, 19))
    record.l_seq = unsafe_load_integer_le(UInt32, buffer, 21)
    record.next_refid = unsafe_load_integer_le(Int32, buffer, 25)
    record.next_pos = unsafe_load_integer_le(Int32, buffer, 29)
    record.tlen = unsafe_load_integer_le(Int32, buffer, 33)
end

function try_read!(x::BAMReader)::Union{Nothing, XAMError}
    # TODO: Must be a memory of previous errors, such that the user can't
    # simply try to read another record if the previous failed.

    # Check if EOF and return that error if so
    buffer = get_nonempty_buffer(x.reader)
    isnothing(buffer) && return XAMErrors.EndOfFile
    
    record = x.record
    vector = record.data

    # Make room for the fixed bytes
    resize!(vector, BAM_FIXED_BYTES)
    remaining = MemoryView(vector)

    # Load in the first BAM_FIXED_BYTES into the record's vector
    while true
        mn = min(length(buffer), length(remaining))
        @assert mn > 0
        slice = @inbounds buffer[1:mn]
        copyto!(remaining, slice)
        consume(x.reader, mn)
        remaining = @inbounds remaining[mn + 1:end]
        if isempty(remaining)
            break
        else
            buffer = get_nonempty_buffer(x.reader)
            isnothing(buffer) && return XAMErrors.TruncatedFile
        end
    end

    unsafe_load_fixed_fields(record, ImmutableMemoryView(vector))

    # Read all the variable length fields
    total_bytes_to_read = record.block_size - 4

    # Ensure the total field length is not shorter than the sum of the
    # variable field names
    seq_bytes = (record.l_seq + UInt32(1)) >> 1
    minimum_bytes = record.l_read_name + record.l_seq + seq_bytes + 4 * record.n_cigar_op
    total_bytes_to_read < minimum_bytes && return XAMErrors.InvalidBlockLength

    # Read all the remaining bytes, overwriting the fixed field bytes,
    # which we have already moved into the integer fields of the struct
    resize!(vector, total_bytes_to_read)
    remaining = MemoryView(vector)

    while !isempty(remaining)
        buffer = get_nonempty_buffer(x.reader)
        isnothing(buffer) && return XAMErrors.TruncatedFile
        mn = min(length(buffer), length(remaining))
        slice = @inbounds buffer[1:mn]
        copyto!(remaining, slice)
        consume(x.reader, mn)
        remaining = @inbounds remaining[mn + 1:end]
    end

    return nothing
end

end # module XAM2
