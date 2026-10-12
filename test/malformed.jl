# Tests of reading malformed BAM records. These must never throw, and records
# must satisfy their invariants even after a failed read.

const MALFORMED_REFS = [("chr1", 1000), ("chr2", 500)]

# CIGAR op codes in BAM
const OP_CODES = Dict('M' => 0, 'I' => 1, 'D' => 2, 'N' => 3, 'S' => 4, 'H' => 5, 'P' => 6, '=' => 7, 'X' => 8)
encode_cigar(ops) = UInt32[(UInt32(len) << 4) | OP_CODES[op] for (len, op) in ops]

le_bytes(x) = collect(reinterpret(UInt8, [htol(x)]))

# Encode a BAM record with arbitrary, possibly invalid fields.
function raw_record(;
        name::Vector{UInt8} = b"r\0" |> collect,
        cigar::Vector{UInt32} = UInt32[],
        seq::Vector{UInt8} = UInt8[], # 4-bit encoded
        qual::Vector{UInt8} = UInt8[],
        aux::Vector{UInt8} = UInt8[],
        refID::Int32 = Int32(-1),
        pos::Int32 = Int32(-1),
        l_read_name::UInt8 = UInt8(length(name)),
        n_cigar_op::UInt16 = UInt16(length(cigar)),
        flag::UInt16 = UInt16(4),
        l_seq::UInt32 = UInt32(length(qual)),
        next_refID::Int32 = Int32(-1),
        next_pos::Int32 = Int32(-1),
        data::Vector{UInt8} = vcat(name, reinterpret(UInt8, htol.(cigar)), seq, qual, aux),
        block_size::UInt32 = UInt32(32 + length(data)),
    )
    io = IOBuffer()
    for x in (block_size, refID, pos, l_read_name, 0xff, UInt16(4680), n_cigar_op, flag, l_seq, next_refID, next_pos, Int32(0))
        write(io, le_bytes(x))
    end
    write(io, data)
    return take!(io)
end

# The uncompressed BAM header, with the given text and references
function raw_bam_header(; text::String = "", refs = MALFORMED_REFS)
    io = IOBuffer()
    write(io, b"BAM\1", le_bytes(UInt32(ncodeunits(text))), text, le_bytes(UInt32(length(refs))))
    for (name, len) in refs
        write(io, le_bytes(UInt32(ncodeunits(name) + 1)), name, 0x00, le_bytes(UInt32(len)))
    end
    return take!(io)
end

# Write a BAM file with the given raw records, and return its path
function malformed_bam(dir::String, records::Vector{UInt8}...; text::String = "", refs = MALFORMED_REFS)
    path = tempname(dir)
    BGZFWriter(open(path, "w")) do io
        write(io, raw_bam_header(; text, refs))
        for record in records
            write(io, record)
        end
    end
    return path
end

# Compress `data` as a BGZF stream, and return the compressed bytes.
# Concatenating such streams gives a BGZF file with a block boundary between them.
function bgzf_compress(dir::String, data::Vector{UInt8})
    path = tempname(dir)
    BGZFWriter(io -> write(io, data), open(path, "w"); n_workers = 1)
    return read(path)
end

# Corrupt the CRC32 of the first BGZF block in `bytes`
function corrupt_crc!(bytes::Vector{UInt8})
    # BSIZE, the block size minus one, is at zero-based offset 16
    block_size = (Int(bytes[17]) | (Int(bytes[18]) << 8)) + 1
    # The CRC32 is followed by the 4-byte ISIZE at the end of the block
    bytes[block_size - 7] ⊻= 0xff
    return bytes
end

# An IO whose reads fail, like a disk error
mutable struct FailingIO <: IO
    open::Bool
end
Base.isopen(io::FailingIO) = io.open
Base.close(io::FailingIO) = (io.open = false; nothing)
Base.eof(::FailingIO) = error("simulated IO failure")

# Read the first record of a BAM file with the given raw records
function read_first(dir::String, records::Vector{UInt8}...)
    return BAMReader(open(malformed_bam(dir, records...)); n_workers = 1, copy = false) do reader
        XAM2.try_read!(reader)
    end
end

# Check that the record is the empty, unmapped record
function is_reset(record::BAMRecord)
    return query_name(record) === nothing &&
        reference_name(record) === nothing &&
        next_reference_name(record) === nothing &&
        left_position(record) === nothing &&
        XAM2.cigar(record) === nothing &&
        sequence(record) === nothing &&
        isempty(aux_data(record))
end

@testset "BAM header" begin
    mktempdir() do dir
        function open_bam(records::Vector{UInt8}...; kwargs...)
            return BAMReader(open(malformed_bam(dir, records...; kwargs...)); n_workers = 1)
        end
        sq_lines = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"
        record = raw_record(name = collect(b"rec\0"), refID = Int32(1), pos = Int32(0))

        # Text header agreeing with the binary references
        reader = open_bam(; text = "@HD\tVN:1.6\n" * sq_lines)
        @test [sq.sn for sq in get_header(reader).seqs] == ["chr1", "chr2"]
        close(reader)

        # The text may be padded with null bytes, which are not part of the header
        reader = open_bam(; text = sq_lines * "\0\0\0")
        @test length(get_header(reader).seqs) == 2
        @test String(MemoryView(get_header(reader))) == sq_lines
        close(reader)

        # The text need not contain @SQ lines
        reader = open_bam(; text = "@HD\tVN:1.6\n")
        @test isempty(get_header(reader).seqs)
        @test references(reader) == References(["chr1", "chr2"], [1000, 500])
        close(reader)

        # A malformed text header, or @SQ lines disagreeing with the binary references,
        # are returned by `get_header`, but records can still be read, since they use
        # the binary references
        for (text, err) in [
                "@SQ\tSN:chr1\tLN:1000\n" => XAMError(XAMErrors.header_reference_mismatch),
                "@SQ\tSN:chr2\tLN:500\n@SQ\tSN:chr1\tLN:1000\n" => XAMError(XAMErrors.header_reference_mismatch),
                "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:501\n" => XAMError(XAMErrors.header_reference_mismatch),
                sq_lines * "@SQ\tSN:chr3\tLN:1\n" => XAMError(XAMErrors.header_reference_mismatch),
                # The location of errors is the line and column in the text header
                "@HD\tVN:1.6\n@XX\tA:1\n" => XAMError(XAMErrors.header_unknown_record_type; line = 2, column = 2),
                "@SQ\tSN:chr1\n" => XAMError(XAMErrors.header_missing_ln; line = 1, column = 1),
            ]
            reader = open_bam(record; text)
            @test get_header(reader) == err
            @test references(reader) == References(["chr1", "chr2"], [1000, 500])
            (rec,) = collect(reader)
            @test (query_name(rec), reference_name(rec)) == ("rec", "chr2")
            close(reader)
        end

        @testset "try_init" begin
            reader = XAM2.try_init(BAMReader, open(malformed_bam(dir; text = sq_lines)); n_workers = 1)
            @test reader isa BAMReader
            close(reader)

            # A malformed text header does not make construction fail
            reader = XAM2.try_init(BAMReader, open(malformed_bam(dir; text = "@XX")); n_workers = 1)
            @test get_header(reader) isa XAMError
            close(reader)

            # ... unless `strict_header` is set
            io = open(malformed_bam(dir; text = "@XX"))
            @test XAM2.try_init(BAMReader, io; n_workers = 1, strict_header = true) isa XAMError
            @test !isopen(io)
            io = open(malformed_bam(dir; text = "@SQ\tSN:chr1\tLN:1000\n"))
            @test XAM2.try_init(BAMReader, io; n_workers = 1, strict_header = true) ==
                XAMError(XAMErrors.header_reference_mismatch)
            @test_throws XAMError BAMReader(open(malformed_bam(dir; text = "@XX")); n_workers = 1, strict_header = true)
            reader = BAMReader(open(malformed_bam(dir; text = sq_lines)); n_workers = 1, strict_header = true)
            @test get_header(reader) isa SAMHeader
            close(reader)

            # On failure, the error is returned, and the input is closed
            io = open(malformed_bam(dir; refs = [("a", 1), ("a", 2)]))
            @test XAM2.try_init(BAMReader, io; n_workers = 1) == XAMError(XAMErrors.duplicate_reference)
            @test !isopen(io)
            # Also when passing a BGZFReader
            path = tempname(dir)
            BGZFWriter(io -> write(io, b"BAM\1", le_bytes(UInt32(100))), open(path, "w"))
            bgzf = BGZFReader(open(path); n_workers = 1)
            @test XAM2.try_init(BAMReader, bgzf) == XAMError(XAMErrors.header_truncated)
            @test !isopen(bgzf)
            # Malformed BGZF data is returned as an error
            io = IOBuffer(b"not a BAM file")
            err = XAM2.try_init(BAMReader, io; n_workers = 1)
            @test err.kind == XAMErrors.bgzf_error
            @test err.bgzf_error isa BGZFError
            @test !isopen(io)
            path = tempname(dir)
            write(path, corrupt_crc!(bgzf_compress(dir, raw_bam_header(; text = "@HD\tVN:1.6\n"))))
            io = open(path)
            @test error_kind(XAM2.try_init(BAMReader, io; n_workers = 1)) == XAMErrors.bgzf_error
            @test !isopen(io)

            # IO errors are thrown, but the input is still closed
            io = FailingIO(true)
            @test_throws ErrorException("simulated IO failure") XAM2.try_init(BAMReader, io; n_workers = 1)
            @test !isopen(io)
        end
    end
end

@testset "Malformed records" begin
    mktempdir() do dir
        # Virtual offset of the first record. All files from `malformed_bam` have the same
        # header, and the header and first records are in the first BGZF block.
        first_offset = BAMReader(r -> r.offset_first_record, open(malformed_bam(dir)); n_workers = 1)
        @test first_offset.file_offset == 0

        @testset "Valid raw record" begin
            record = read_first(dir, raw_record(refID = Int32(1), pos = Int32(9), next_refID = Int32(0), next_pos = Int32(0)))
            @test record isa BAMRecord
            @test query_name(record) == "r"
            @test reference_name(record) == "chr2"
            @test next_reference_name(record) == "chr1"
            @test left_position(record) == 10
            @test quality_scores(record) === nothing

            # Undefined flag bits are cleared
            record = read_first(dir, raw_record(flag = 0xffff))
            @test flags(record) == SAMFlags(0xffff)
            @test Integer(flags(record)) == 0x0fff
        end

        @testset "Structural errors" begin
            for (kwargs, err) in [
                    (; block_size = UInt32(31)) => XAMErrors.invalid_block_length,
                    (; block_size = UInt32(0)) => XAMErrors.invalid_block_length,
                    # Fields don't fit in block
                    (; l_read_name = UInt8(3)) => XAMErrors.invalid_block_length,
                    (; n_cigar_op = UInt16(1)) => XAMErrors.invalid_block_length,
                    (; l_seq = UInt32(1)) => XAMErrors.invalid_block_length,
                    # Regression: These overflowed when computed as UInt32
                    (; l_seq = typemax(UInt32)) => XAMErrors.invalid_block_length,
                    (; l_seq = 0xaaaaaaab) => XAMErrors.invalid_block_length,
                    # Read name
                    (; l_read_name = UInt8(0)) => XAMErrors.invalid_read_name,
                    (; name = collect(b"ab")) => XAMErrors.invalid_read_name,
                    # Reference IDs
                    (; refID = Int32(2)) => XAMErrors.invalid_reference_id,
                    (; refID = Int32(-2)) => XAMErrors.invalid_reference_id,
                    (; next_refID = Int32(2)) => XAMErrors.invalid_reference_id,
                    (; next_refID = typemin(Int32)) => XAMErrors.invalid_reference_id,
                    # Positions
                    (; pos = Int32(-2)) => XAMErrors.invalid_position,
                    (; next_pos = typemin(Int32)) => XAMErrors.invalid_next_position,
                    # SAM can't represent the 1-based position 2^31
                    (; pos = typemax(Int32)) => XAMErrors.invalid_position,
                    (; next_pos = typemax(Int32)) => XAMErrors.invalid_next_position,
                ]
                @test read_first(dir, raw_record(; kwargs...)) == XAMError(err; offset = first_offset)
            end

            # The largest valid position, which fits in SAM
            record = read_first(dir, raw_record(refID = Int32(0), pos = typemax(Int32) - Int32(1)))
            @test left_position(record) == typemax(Int32)
            @test left_position(XAM2.try_copy!(SAMRecord(), record)) == typemax(Int32)

            # Truncated in fixed fields, and in variable-length fields
            @test read_first(dir, raw_record()[1:20]) == XAMError(XAMErrors.truncated_file; offset = first_offset)
            @test read_first(dir, raw_record(; qual = fill(0x1e, 10), seq = fill(0x11, 5))[1:40]) ==
                XAMError(XAMErrors.truncated_file; offset = first_offset)
        end

        @testset "Errors are sticky, and the record stays valid" begin
            good = raw_record(name = collect(b"good\0"), refID = Int32(0), pos = Int32(0))
            bad = raw_record(refID = Int32(5), pos = Int32(0), name = collect(b"bad\0"))
            path = malformed_bam(dir, good, bad)
            bad_offset = VirtualOffset(0, first_offset.block_offset + length(good))
            expected = XAMError(XAMErrors.invalid_reference_id; offset = bad_offset)
            @test sprint(showerror, expected) ==
                "XAMError: record's reference ID is not in the references at virtual offset $(bad_offset)"
            BAMReader(open(path); n_workers = 1, copy = false) do reader
                record = XAM2.try_read!(reader)
                @test query_name(record) == "good"
                @test XAM2.try_read!(reader) == expected
                # The record previously returned is reset, not partially overwritten
                @test is_reset(record)
                @test XAM2.try_read!(reader) == expected
                @test_throws XAMError iterate(reader)

                # Seeking clears the error
                seek_first_record(reader)
                @test query_name(XAM2.try_read!(reader)) == "good"
                @test XAM2.try_read!(reader) == expected

                # Seeking to the offset of the error reproduces it
                seek(reader, first_offset)
                @test query_name(XAM2.try_read!(reader)) == "good"
                seek(reader, expected.offset)
                @test XAM2.try_read!(reader) == expected
            end
        end

        @testset "Malformed BGZF" begin
            header = raw_bam_header()
            rec1 = raw_record(name = collect(b"rec1\0"), refID = Int32(0), pos = Int32(0))
            rec2 = raw_record(name = collect(b"rec2\0"), refID = Int32(1), pos = Int32(0))

            # Corrupt block at a record boundary: The error's offset is the start of the
            # malformed block, where the next record starts.
            first_stream = bgzf_compress(dir, vcat(header, rec1))
            path = tempname(dir)
            write(path, vcat(first_stream, corrupt_crc!(bgzf_compress(dir, rec2))))
            BAMReader(open(path); n_workers = 1, copy = false) do reader
                @test query_name(XAM2.try_read!(reader)) == "rec1"
                err = XAM2.try_read!(reader)
                @test err.kind == XAMErrors.bgzf_error
                @test err.bgzf_error isa BGZFError
                @test err.bgzf_error.file_offset == length(first_stream)
                @test err.offset == VirtualOffset(length(first_stream), 0)
                @test startswith(sprint(showerror, err), "XAMError: malformed BGZF compression (BGZFError: ")
                # Sticky
                @test XAM2.try_read!(reader) === err
                @test_throws XAMError iterate(reader)
                # Seeking to the offset reproduces the error, without throwing
                seek(reader, err.offset)
                @test XAM2.try_read!(reader) == err
                # Seeking back recovers
                seek_first_record(reader)
                @test query_name(XAM2.try_read!(reader)) == "rec1"
                @test XAM2.try_read!(reader) == err
            end

            # Corrupt block in the middle of a record: The offset is the record's start
            path = tempname(dir)
            first_stream = bgzf_compress(dir, vcat(header, rec1, rec2[1:10]))
            write(path, vcat(first_stream, corrupt_crc!(bgzf_compress(dir, rec2[11:end]))))
            BAMReader(open(path); n_workers = 1, copy = false) do reader
                @test query_name(XAM2.try_read!(reader)) == "rec1"
                err = XAM2.try_read!(reader)
                @test err.kind == XAMErrors.bgzf_error
                @test err.offset == VirtualOffset(0, length(header) + length(rec1))
                @test is_reset(reader.record)
                seek(reader, err.offset)
                @test XAM2.try_read!(reader) == err
            end

            # Seeking outside a block is invalid input, not malformed data, so it throws
            BAMReader(open(malformed_bam(dir, rec1)); n_workers = 1) do reader
                @test_throws BGZFError seek(reader, VirtualOffset(0, 60000))
            end
        end

        @testset "Malformed content" begin
            # Invalid CIGAR op: Accessors return the error, and don't throw
            record = read_first(dir, raw_record(refID = Int32(0), pos = Int32(0), cigar = UInt32[(4 << 4) | 9], seq = [0x12, 0x48], qual = fill(0x1e, 4)))
            @test record isa BAMRecord
            @test XAM2.cigar(record) isa CIGARError
            @test mapping_span(record) isa CIGARError
            @test sprint(show, MIME"text/plain"(), record) isa String

            # CIGAR query length must match the sequence length
            seq_fields = (; refID = Int32(0), pos = Int32(0), seq = [0x12, 0x48], qual = fill(0x1e, 4))
            record = read_first(dir, raw_record(; seq_fields..., cigar = encode_cigar([(3, 'M')])))
            @test XAM2.cigar(record) == XAMError(XAMErrors.cigar_query_length_mismatch)
            @test mapping_span(record) == XAMError(XAMErrors.cigar_query_length_mismatch)
            @test sprint(show, MIME"text/plain"(), record) isa String
            # Hard clips are not part of the query
            record = read_first(dir, raw_record(; seq_fields..., cigar = encode_cigar([(2, 'H'), (4, 'M'), (1, 'H')])))
            @test XAM2.cigar(record) isa XAM2.BAMCIGAR
            # Any CIGAR is allowed if the sequence is missing
            record = read_first(dir, raw_record(; refID = Int32(0), pos = Int32(0), cigar = encode_cigar([(10, 'M')])))
            @test XAM2.cigar(record) isa XAM2.BAMCIGAR
            @test mapping_span(record) == 1:10

            # Missing sequence: Quality is also missing
            record = read_first(dir, raw_record(refID = Int32(0), pos = Int32(0)))
            @test sequence(record) === nothing
            @test quality_scores(record) === nothing
        end

        @testset "CIGAR in CG tag" begin
            placeholder = encode_cigar([(4, 'S'), (1000, 'N')])
            real_cigar = encode_cigar([(2, 'M'), (1, 'I'), (1, 'M')])
            cg_tag = vcat(collect(b"CGBI"), le_bytes(UInt32(3)), reinterpret(UInt8, htol.(real_cigar)))
            fields = (; refID = Int32(0), pos = Int32(0), cigar = placeholder, seq = [0x12, 0x48], qual = fill(0x1e, 4))

            # The CIGAR is loaded from the CG tag
            record = read_first(dir, raw_record(; fields..., aux = cg_tag))
            @test String(MemoryView(CIGAR(XAM2.cigar(record)))) == "2M1I1M"
            @test mapping_span(record) == 1:3

            # Like htslib writes it, m in the placeholder `kSmN` is the reference length
            # of the alignment, not of the reference sequence
            alignment_placeholder = encode_cigar([(4, 'S'), (3, 'N')])
            record = read_first(dir, raw_record(; fields..., cigar = alignment_placeholder, aux = cg_tag))
            @test String(MemoryView(CIGAR(XAM2.cigar(record)))) == "2M1I1M"

            # The CIGAR in the CG tag must also match the sequence length
            bad_cigar = encode_cigar([(2, 'M'), (1, 'I'), (2, 'M')])
            bad_tag = vcat(collect(b"CGBI"), le_bytes(UInt32(3)), reinterpret(UInt8, htol.(bad_cigar)))
            record = read_first(dir, raw_record(; fields..., aux = bad_tag))
            @test XAM2.cigar(record) == XAMError(XAMErrors.cigar_query_length_mismatch)

            # Without a CG tag, or with malformed aux data, the placeholder is used
            for aux in [UInt8[], collect(b"CGB"), collect(b"XYZabc")]
                record = read_first(dir, raw_record(; fields..., aux))
                @test String(MemoryView(CIGAR(XAM2.cigar(record)))) == "4S1000N"
            end

            # Unmapped records can't have the placeholder CIGAR, since they have no reference.
            # Regression: This threw a BoundsError.
            record = read_first(dir, raw_record(; fields..., refID = Int32(-1), aux = cg_tag))
            @test String(MemoryView(CIGAR(XAM2.cigar(record)))) == "4S1000N"
        end
    end
end
