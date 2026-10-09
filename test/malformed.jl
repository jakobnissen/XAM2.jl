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

# Write a BAM file with the given raw records, and return its path
function malformed_bam(dir::String, records::Vector{UInt8}...)
    path = tempname(dir)
    BGZFWriter(open(path, "w")) do io
        write(io, b"BAM\1", le_bytes(UInt32(0)), le_bytes(UInt32(length(MALFORMED_REFS))))
        for (name, len) in MALFORMED_REFS
            write(io, le_bytes(UInt32(ncodeunits(name) + 1)), name, 0x00, le_bytes(UInt32(len)))
        end
        for record in records
            write(io, record)
        end
    end
    return path
end

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

@testset "Malformed records" begin
    mktempdir() do dir
        @testset "Valid raw record" begin
            record = read_first(dir, raw_record(refID = Int32(1), pos = Int32(9), next_refID = Int32(0), next_pos = Int32(0)))
            @test record isa BAMRecord
            @test query_name(record) == "r"
            @test reference_name(record) == "chr2"
            @test next_reference_name(record) == "chr1"
            @test left_position(record) == 10
            @test quality_scores(record) === nothing
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
                    (; next_pos = typemin(Int32)) => XAMErrors.invalid_position,
                ]
                @test read_first(dir, raw_record(; kwargs...)) == err
            end

            # Truncated in fixed fields, and in variable-length fields
            @test read_first(dir, raw_record()[1:20]) == XAMErrors.truncated_file
            @test read_first(dir, raw_record(; qual = fill(0x1e, 10), seq = fill(0x11, 5))[1:40]) ==
                XAMErrors.truncated_file
        end

        @testset "Errors are sticky, and the record stays valid" begin
            good = raw_record(name = collect(b"good\0"), refID = Int32(0), pos = Int32(0))
            bad = raw_record(refID = Int32(5), pos = Int32(0), name = collect(b"bad\0"))
            path = malformed_bam(dir, good, bad)
            BAMReader(open(path); n_workers = 1, copy = false) do reader
                record = XAM2.try_read!(reader)
                @test query_name(record) == "good"
                @test XAM2.try_read!(reader) == XAMErrors.invalid_reference_id
                # The record previously returned is reset, not partially overwritten
                @test is_reset(record)
                @test XAM2.try_read!(reader) == XAMErrors.invalid_reference_id
                @test_throws XAMError iterate(reader)

                # Seeking clears the error
                seek_first_record(reader)
                @test query_name(XAM2.try_read!(reader)) == "good"
                @test XAM2.try_read!(reader) == XAMErrors.invalid_reference_id
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
            @test XAM2.cigar(record) == XAMErrors.cigar_query_length_mismatch
            @test mapping_span(record) == XAMErrors.cigar_query_length_mismatch
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

            # The CIGAR in the CG tag must also match the sequence length
            bad_cigar = encode_cigar([(2, 'M'), (1, 'I'), (2, 'M')])
            bad_tag = vcat(collect(b"CGBI"), le_bytes(UInt32(3)), reinterpret(UInt8, htol.(bad_cigar)))
            record = read_first(dir, raw_record(; fields..., aux = bad_tag))
            @test XAM2.cigar(record) == XAMErrors.cigar_query_length_mismatch

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
