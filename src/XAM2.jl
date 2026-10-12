module XAM2

using BufferIO: BufferIO,
    AbstractBufReader,
    AbstractBufWriter,
    BufReader,
    VecWriter,
    get_nonempty_buffer,
    get_unflushed,
    consume,
    fill_buffer,
    get_buffer,
    grow_buffer,
    line_views

using MemoryViews:
    MemoryViews,
    ImmutableMemoryView,
    MemoryView,
    MutableMemoryView,
    split_last,
    split_each

using StringViews: StringView
using BGZFLib: BGZFReader, BGZFWriter, BGZFError, BGZFErrors, virtual_position, VirtualOffset
using CIGARStrings: CIGARStrings, AbstractCIGAR, BAMCIGAR, CIGAR, CIGARError, unsafe_switch_memory
using XAMAuxData: XAMAuxData, AuxTag, try_auxtag, is_well_formed

using Automa: Automa, @re_str, onenter!, onexit!, onall!, rep, rep1, opt

export BAMRecord,
    BAMReader,
    BAMWriter,
    SAMRecord,
    SAMReader,
    SAMHeader,
    SAMFlag,
    SAMFlags,
    Flag,
    XAMRecord,
    push,
    StringVec,
    References,
    references,
    get_header,
    query_name,
    flags,
    reference_name,
    reference_length,
    next_reference_name,
    reference_index,
    next_reference_index,
    left_position,
    next_position,
    mapping_quality,
    template_length,
    sequence,
    quality_string,
    quality_scores,
    aux_data,
    is_mapped,
    is_first,
    is_last,
    mapping_span,
    XAMError,
    XAMErrorKind,
    XAMErrors

public try_parse!,
    try_read!,
    cigar,
    seek_first_record,
    try_init,
    try_copy!,
    QualityView,
    try_write,
    Reference,
    reference,
    try_parse_header,
    HasEq

"""
    module XAMErrors

This module is used as a namespace for the [`XAMErrorKind`](@ref) enum.
"""
module XAMErrors

    """
        XAMErrorKind

    Enum that contains information about a [`XAMError`](@ref).
    New variants may be added in future, minor releases.
    The integer value of these kinds are not stable over minor releases.

    julia> kind = XAMErrors.invalid_cigar;

    julia> throw(XAMError(kind))
    ERROR:
    [...]
    """
    @enum XAMErrorKind::UInt8 begin
        # IO
        truncated_file
        closed_stream
        bgzf_error # malformed BGZF compression, see the `bgzf_error` property of `XAMError`
        # BAM header
        header_truncated
        header_bad_magic_bytes
        header_top_bit_set
        header_reference_mismatch # @SQ lines in BAM header differ from binary references
        # SAM header text, also in BAM
        header_invalid_line # line is not `@`, a two-letter record type, and a tab
        header_unknown_record_type
        header_invalid_field # field is not `TAG:VALUE` with a valid tag and non-empty value
        header_duplicate_tag # tag occurs twice in the same line
        header_hd_not_first
        header_missing_vn
        header_invalid_vn
        header_invalid_sort_order
        header_missing_sn
        header_invalid_sn
        header_duplicate_sn
        header_missing_ln
        header_invalid_ln
        header_invalid_topology
        header_missing_id # in @RG or @PG line
        # References
        invalid_reference # reference with invalid name or length
        duplicate_reference # two references with the same name
        invalid_reference_id # record has a reference ID outside the references
        unknown_reference # record's reference is not present in the destination
        reference_length_mismatch # same reference name, but different lengths
        # Record structure
        invalid_block_length
        record_too_long # BAM record's block size, or a CIGAR placeholder length, does not fit
        missing_fields # SAM line has fewer than 11 fields
        line_too_long # more than 4 GiB line
        qual_no_seq # SAM record has a quality, but no sequence
        qual_seq_mismatch # SAM record's quality and sequence have different lengths
        # Record fields
        invalid_read_name # BAM: length is zero, or not null terminated. SAM: invalid QNAME
        invalid_flags # SAM FLAG
        invalid_reference_name # SAM RNAME
        invalid_position # position is less than -1 (BAM), or invalid POS (SAM)
        invalid_mapq # SAM MAPQ
        invalid_cigar
        invalid_next_reference_name # SAM RNEXT
        invalid_next_position # next position is less than -1 (BAM), or invalid PNEXT (SAM)
        invalid_template_length # SAM TLEN
        invalid_sequence # SAM SEQ
        invalid_quality # SAM QUAL
        invalid_aux_data # SAM auxiliary fields
        # Record content, validated lazily
        cigar_query_length_mismatch # CIGAR's query length differs from sequence length
    end
end # module XAMErrors

using .XAMErrors: XAMErrorKind

"""
    XAMError <: Exception

An error from reading, parsing or writing invalid SAM or BAM data.
Functions such as `XAM2.try_read!` return an `XAMError`, and their throwing
counterparts, such as iteration, throw it.

An `XAMError` has the following public properties:
* `kind::XAMErrorKind`: What went wrong. Use this to check the kind of error
  programmatically, see [`XAMErrors`](@ref).
* `record::Union{Nothing, Int}`: The 1-based number of the record in the file
  where the error occurred, counted from the first record.
* `line::Union{Nothing, Int}`: The 1-based line number, for errors in SAM text,
  i.e. SAM files, and the text headers of SAM and BAM files.
* `column::Union{Nothing, Int}`: The 1-based byte column in the line.
* `offset::Union{Nothing, VirtualOffset}`: For BAM records, the BGZF virtual offset
  of the start of the record. Seeking a `BAMReader` to this offset and reading
  reproduces the error.
* `bgzf_error::Union{Nothing, BGZFLib.BGZFError}`: For errors of kind
  `XAMErrors.bgzf_error`, the error from decompressing the malformed BGZF data.

The location properties are `nothing` if not applicable or not known.
SAM records have `record`, `line` and `column`; text headers of SAM and BAM files have
`line` and `column`; and BAM records have `offset`.

    XAMError(
        kind::XAMErrorKind;
        record=nothing, line=nothing, column=nothing, offset=nothing, bgzf_error=nothing,
    )
"""
struct XAMError <: Exception
    kind::XAMErrorKind
    record::Union{Nothing, Int}
    line::Union{Nothing, Int}
    column::Union{Nothing, Int}
    offset::Union{Nothing, VirtualOffset}
    bgzf_error::Union{Nothing, BGZFError}
end

function XAMError(
        kind::XAMErrorKind;
        record::Union{Nothing, Int} = nothing,
        line::Union{Nothing, Int} = nothing,
        column::Union{Nothing, Int} = nothing,
        offset::Union{Nothing, VirtualOffset} = nothing,
        bgzf_error::Union{Nothing, BGZFError} = nothing,
    )
    return XAMError(kind, record, line, column, offset, bgzf_error)
end

XAMError(e::BGZFError; offset::Union{Nothing, VirtualOffset} = nothing) = XAMError(XAMErrors.bgzf_error; offset, bgzf_error = e)


function error_message(kind::XAMErrorKind)::String
    E = XAMErrors
    return if kind == E.truncated_file
        "input data is truncated; does not contain whole record or header"
    elseif kind == E.closed_stream
        "stream is closed"
    elseif kind == E.bgzf_error
        "malformed BGZF compression"
    elseif kind == E.header_truncated
        "BAM header is truncated"
    elseif kind == E.header_bad_magic_bytes
        "data does not begin with the BAM magic bytes"
    elseif kind == E.header_top_bit_set
        "header text length or number of references in BAM header exceeds 2147483647"
    elseif kind == E.header_reference_mismatch
        "@SQ lines in BAM header do not match the reference sequences in the binary header"
    elseif kind == E.header_invalid_line
        "SAM header line does not start with '@' followed by a two-letter record type and a tab"
    elseif kind == E.header_unknown_record_type
        "SAM header record type is not one of HD, SQ, RG, PG or CO"
    elseif kind == E.header_invalid_field
        "SAM header field is not of the form TAG:VALUE, with a tag matching " *
            "[A-Za-z][A-Za-z0-9] and a non-empty value"
    elseif kind == E.header_duplicate_tag
        "tag occurs more than once in the same SAM header line"
    elseif kind == E.header_hd_not_first
        "@HD line is not the first line of the SAM header"
    elseif kind == E.header_missing_vn
        "@HD line has no VN field"
    elseif kind == E.header_invalid_vn
        "VN field of @HD line is not of the form MAJOR.MINOR"
    elseif kind == E.header_invalid_sort_order
        "SO field of @HD line is not \"unknown\", \"unsorted\", \"queryname\" or \"coordinate\""
    elseif kind == E.header_missing_sn
        "@SQ line has no SN field"
    elseif kind == E.header_invalid_sn
        "SN field of @SQ line is not a valid reference name"
    elseif kind == E.header_duplicate_sn
        "SN field of @SQ line is the same as that of a previous @SQ line"
    elseif kind == E.header_missing_ln
        "@SQ line has no LN field"
    elseif kind == E.header_invalid_ln
        "LN field of @SQ line is not an integer in 1:2147483647"
    elseif kind == E.header_invalid_topology
        "TP field of @SQ line is not \"linear\" or \"circular\""
    elseif kind == E.header_missing_id
        "@RG or @PG line has no ID field"
    elseif kind == E.invalid_reference
        "reference sequence has an invalid name or length"
    elseif kind == E.duplicate_reference
        "two reference sequences have the same name"
    elseif kind == E.invalid_reference_id
        "record's reference ID is not in the references"
    elseif kind == E.unknown_reference
        "record's reference is not present in the destination's references"
    elseif kind == E.reference_length_mismatch
        "record's reference has a different length in the destination's references"
    elseif kind == E.invalid_block_length
        "BAM record's block size is too small to contain its fields"
    elseif kind == E.record_too_long
        "record is too long to be represented in BAM"
    elseif kind == E.missing_fields
        "SAM record has fewer than 11 tab-separated fields"
    elseif kind == E.line_too_long
        "SAM line is longer than 4294967294 bytes"
    elseif kind == E.qual_no_seq
        "SAM record has a quality (QUAL), but no sequence (SEQ)"
    elseif kind == E.qual_seq_mismatch
        "SAM record's quality (QUAL) and sequence (SEQ) have different lengths"
    elseif kind == E.invalid_read_name
        "invalid read name (QNAME)"
    elseif kind == E.invalid_flags
        "invalid flags (FLAG)"
    elseif kind == E.invalid_reference_name
        "invalid reference name (RNAME)"
    elseif kind == E.invalid_position
        "invalid position (POS)"
    elseif kind == E.invalid_mapq
        "invalid mapping quality (MAPQ)"
    elseif kind == E.invalid_cigar
        "invalid CIGAR"
    elseif kind == E.invalid_next_reference_name
        "invalid next reference name (RNEXT)"
    elseif kind == E.invalid_next_position
        "invalid next position (PNEXT)"
    elseif kind == E.invalid_template_length
        "invalid template length (TLEN)"
    elseif kind == E.invalid_sequence
        "invalid sequence (SEQ)"
    elseif kind == E.invalid_quality
        "invalid quality (QUAL)"
    elseif kind == E.invalid_aux_data
        "invalid auxiliary data"
    else
        @assert kind == E.cigar_query_length_mismatch
        "query length of CIGAR differs from the length of the sequence"
    end
end

function Base.showerror(io::IO, error::XAMError)
    print(io, "XAMError: ", error_message(error.kind))
    bgzf_error = error.bgzf_error
    if bgzf_error !== nothing
        print(io, " (")
        showerror(io, bgzf_error)
        print(io, ')')
    end
    is_first = true
    for (name, value) in (("record", error.record), ("line", error.line), ("column", error.column))
        value === nothing && continue
        print(io, is_first ? " at " : ", ", name, ' ', value)
        is_first = false
    end
    offset = error.offset
    if offset !== nothing
        print(io, is_first ? " at " : ", ", "virtual offset ", offset)
    end
    return nothing
end

struct Uninit end

"""
    HasEq()

Zero-sized struct indicating that a BAM sequence contains `=`, and therefore cannot
be returned as a `LongSequence`.
"""
struct HasEq end

const CTX = Automa.CodeGenContext(; generator = :goto, clean = true)

@static if isdefined(Base, :memoryindex)
    memory_index(x::MemoryRef) = Base.memoryindex(x)
else
    memory_index(x::MemoryRef) = Core.memoryrefoffset(x)
end

# Index in underlying memory of first index of MemoryView or Vector
memory_index(v::Vector) = memory_index(Base.cconvert(Ptr, v))
memory_index(m::MemoryView) = first(only(parentindices(m)))

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

function tryparse_base10(::Type{T}, mem::ImmutableMemoryView{UInt8}) where {T <: Union{Int8, Int16, Int32}}
    isempty(mem) && return nothing
    start = @inbounds mem[1]
    (is_positive, maxval, mem) = if start == UInt8('-')
        length(mem) < 2 && return nothing
        (false, abs(Int64(typemin(T))) % UInt64, @inbounds(mem[2:end]))
    else
        (true, typemax(T) % UInt64, mem)
    end
    n = @inline tryparse_base10(UInt64, mem)
    isnothing(n) && return nothing
    n > maxval && return nothing
    n = n % T
    return is_positive ? n : -n
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

function inbounds_write_base10!(mem::MutableMemoryView{UInt8}, n::Unsigned)::Int
    T = typeof(n)
    if iszero(n)
        @inbounds mem[1] = 0x30
        return 1
    end
    n_digits = 0
    @inbounds while !iszero(n)
        n_digits += 1
        (n, rem) = divrem(n, T(10))
        mem[n_digits] = (rem % UInt8) + 0x30
    end
    rv = n_digits
    @inbounds for i in 1:(n_digits >> 1)
        (mem[i], mem[rv]) = (mem[rv], mem[i])
        rv -= 1
    end
    return n_digits
end

# The BAI bin of the 0-based, half-open region `beg:(stop - 1)`, see the SAM specification,
# section 5.3. The region must not be empty, i.e. `stop > beg`. Unmapped records with no
# position, which the specification bins as `reg2bin(-1, 0)`, have bin 4680.
# Positions beyond the range of BAI (2^29) give bins that don't fit in 16 bits,
# which are truncated, like htslib does.
function reg2bin(beg::T, stop::T)::UInt16 where {T <: Unsigned}
    stop -= one(T)
    # Truncating the terms before adding them gives the same result as truncating the sum
    return if beg >> 14 == stop >> 14
        (div((1 << 15) - 1, 7) % UInt16) + ((beg >> 14) % UInt16)
    elseif beg >> 17 == stop >> 17
        (div((1 << 12) - 1, 7) % UInt16) + ((beg >> 17) % UInt16)
    elseif beg >> 20 == stop >> 20
        (div((1 << 9) - 1, 7) % UInt16) + ((beg >> 20) % UInt16)
    elseif beg >> 23 == stop >> 23
        (div((1 << 6) - 1, 7) % UInt16) + ((beg >> 23) % UInt16)
    elseif beg >> 26 == stop >> 26
        (div((1 << 3) - 1, 7) % UInt16) + ((beg >> 26) % UInt16)
    else
        0x0000
    end
end

function reg2bins!(bins::Vector{UInt16}, span::UnitRange{Int32})
    (start, stop) = (first(span) - Int32(1), last(span) - Int32(1))
    # @assert start in Int32(1):Int32(2^29 - 1)
    # @assert stop in start:(Int32(2^29 - 1))

    # If start is zero, the read is unmapped, and so the range
    # must be 0:0 exactly
    # @assert !iszero(start) || iszero(stop)

    # TODO: Check the 0,0 case
    idxs = (1, 9, 73, 585, 4681)
    shifts = (26, 23, 20, 17, 14)

    ranges = map(1:5) do i
        idx = @inbounds(idxs[i]) % UInt16
        shift = @inbounds shifts[i] & 63
        (idx + (start >> shift) % UInt16):(idx + (stop >> shift) % UInt16)
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

# This internal function _assumes_ that:
# * If old_cigar is a AbstractCIGAR, it aliases `old_vector`
# * old_vector == new_vector
function copy_cigar!(
        old_cigar::Union{Nothing, Uninit, CIGARError, XAMError, AbstractCIGAR},
        old_vector::Vector{UInt8},
        new_vector::Vector{UInt8},
    )
    if !isa(old_cigar, AbstractCIGAR)
        return old_cigar
    end
    # First: Which index in `old_vector` does the backing memory of
    # `old_cigar` start?
    cigar_mem = MemoryView(old_cigar)
    index_cigar_mem = memory_index(cigar_mem)
    index_old_vec = memory_index(old_vector)
    cigar_mem_start = index_cigar_mem - index_old_vec + 1

    # Now, copy the bytes of cigar_mem into `new_vector` as those same
    # indices
    # TODO: Inbounds once we've tested correctness
    dst_mem = MemoryView(new_vector)[cigar_mem_start:((cigar_mem_start + length(cigar_mem)) - 1)]
    copy!(dst_mem, cigar_mem) # TODO: Inbounds once we've checked
    dst_mem = ImmutableMemoryView(dst_mem)
    new_cigar = unsafe_switch_memory(old_cigar, dst_mem)
    return new_cigar
end

# Return `cig`, or an error if its query length differs from the sequence length.
# A sequence length of zero means the sequence is missing, so any CIGAR is allowed.
function check_query_length(cig::T, seq_len::Int)::Union{T, XAMError} where {T <: AbstractCIGAR}
    # The query length is stored in the CIGAR, so this is cheap
    if !iszero(seq_len) && CIGARStrings.query_length(cig) != seq_len
        return XAMError(XAMErrors.cigar_query_length_mismatch)
    end
    return cig
end

"""
    seek_first_record(reader::Union{SAMReader, BAMReader}) -> reader

Seek `reader` back to the first record.
"""
function seek_first_record end

# Run `f()`, closing `io` if it throws. Used when constructing readers, which take
# ownership of `io`: If construction fails, the caller has no reader to close.
# The exception is rethrown, not handled.
function close_on_throw(f, io)
    return try
        f()
    catch
        close(io)
        rethrow()
    end
end

# Errors are sticky until the reader is seeked, except if it's closed
function clear_error!(reader)
    err = reader.error
    if err !== nothing && err.kind != XAMErrors.closed_stream
        reader.error = nothing
    end
    return nothing
end

include("flags.jl")
include("stringvec.jl")
include("header.jl")
include("references.jl")
include("bamseq.jl")
include("qualityview.jl")
include("sam.jl")
include("bam.jl")

"""
Alias for Union{SAMRecord, BAMRecord}
"""
const XAMRecord = Union{SAMRecord, BAMRecord}

"""
Alias for Union{SAMReader, BAMReader}
"""
const XAMReader = Union{SAMReader, BAMReader}

"""
    query_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the query name (QNAME) of `record`, or `nothing` if the data
is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

!!! warning
    The returned string aliases the record. Mutating the record will
    invalidate the string.
"""
function query_name end

"""
    reference_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the reference name (RNAME) of `record`, or `nothing` if the data
is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

For `SAMRecord`, this function does _not_ validate that the reference
name is present in the SAM header.

!!! warning
    For `SAMRecord`, this aliases the record, whereas for `BAMRecord`,
    the data aliases the `References` of the record.
    Mutating the underlying record will invalidate the string if returned
    from a `SAMRecord`
"""
function reference_name end

"""
    next_reference_name(record::XAMRecord)::Union{Nothing, AbstractString}

Return the reference name of the next segment in the template (RNEXT) of `record`,
or `nothing` if the data is not available.

The name is returned as a UTF-8 encoded `AbstractString`.

!!! warning
    For `SAMRecord`, this aliases the record, whereas for `BAMRecord`,
    the data aliases the `References` of the record.
    Mutating the underlying record will invalidate the string if returned
    from a `SAMRecord`
"""
function next_reference_name end

"""
    reference_length(record::BAMRecord)::Union{Nothing, Int}

Return the length of the reference that `record` is mapped to, or `nothing`
if the record is not mapped.
"""
function reference_length end

"""
    left_position(record::XAMRecord)::Union{Nothing, Int}

Return the 1-based leftmost (i.e. 5' of the reference) mapped position
of `record`, or `nothing` if the data is not available.
"""
function left_position end

"""
    next_position(record::XAMRecord)::Union{Nothing, Int}

Return the 1-based leftmost position of the next segment in the template (PNEXT)
of `record`, or `nothing` if the data is not available.

See also: [`next_reference_name`](@ref), [`left_position`](@ref)
"""
function next_position end

"""
    mapping_quality(record::XAMRecord)::Union{Nothing, Int}

Return the mapping quality (MAPQ) of `record`, in `0:254`, or `nothing` if it is
not available (stored as 255).
"""
function mapping_quality end

"""
    template_length(record::XAMRecord)::Union{Nothing, Int}

Return the observed template length (TLEN) of `record`, or `nothing` if it is
not available (stored as 0). It is positive for the leftmost segment of the template,
and negative for the rightmost.
"""
function template_length end

"""
    cigar(record::SAMRecord)::Union{Nothing, CIGAR, CIGARError, XAMError}
    cigar(record::BAMRecord)::Union{Nothing, BAMCIGAR, CIGARError, XAMError}

Get the `cigar` of `record` as a `CIGARStrings.CIGAR` or
`CIGARStrings.BAMCIGAR`, or return `nothing` if the data is not available,
or a `CIGARStrings.CIGARError` if the cigar is malformed.
If the record has a sequence, and the query length of the cigar differs from
the length of the sequence, return an `XAMError` of kind
`XAMErrors.cigar_query_length_mismatch`.

This operation will validate the underlying data, and then cache
the constructed cigar such that subsequent calls to this function
is almost free.

!!! warning
    The returned `AbstractCIGAR` object aliases the record.
    Mutating the record will invalidate the cigar.
"""
function cigar end

"""
    flags(record::XAMRecord)::SAMFlags

Get the flags of `record`.
"""
flags(x::XAMRecord) = x.flags

"""
    is_mapped(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.UNMAP`.
"""
is_mapped(x::XAMRecord) = !in(UNMAP, flags(x))

"""
    is_first(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.READ1`.
"""
is_first(x::XAMRecord) = in(READ1, flags(x))

"""
    is_last(record::XAMRecord)::Bool

Check if the flags of `record` contain `SAMFlag.READ2`.
"""
is_last(x::XAMRecord) = in(READ2, flags(x))

"""
    get_header(::SAMReader)::SAMHeader
    get_header(::BAMReader)::Union{SAMHeader, XAMError}

Get the header of the reader.

For a `BAMReader`, the reference sequences are read from the binary part of the
header (see [`references`](@ref)), so the text header is not needed to read records.
If the text header is malformed, or its `@SQ` lines do not match the binary references,
this returns an `XAMError`, and the records can still be read.
"""
function get_header end

"""
    mapping_span(record::XAMRecord)::Union{Nothing, UnitRange{Int}, CIGARError, XAMError}

Return the 1-based positions in the reference that `record` is mapped to,
or `nothing` if the position or CIGAR is not available, or the error
returned by [`cigar`](@ref) if the CIGAR is invalid.

E.g. a record with CIGAR `10M100S` mapping to reference position 5
would map to `5:14`.
"""
function mapping_span(record::XAMRecord)::Union{Nothing, UnitRange{Int}, CIGARError, XAMError}
    left = @something left_position(record) return nothing
    cig = cigar(record)
    cig === nothing && return nothing
    cig isa Union{CIGARError, XAMError} && return cig
    return left:(left + CIGARStrings.ref_length(cig) - 1)
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

end # module XAM2
