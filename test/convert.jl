# Tests of converting between SAMRecord and BAMRecord

# The fields of a BAM record that are stored in the BAM file. Auxiliary values are
# compared by value, since integers may be encoded with different types.
function bam_fields(r::BAMRecord)
    fixed = Tuple(getfield(r, f) for f in fieldnames(BAMRecord) if f ∉ (:cigar, :refs, :data))
    return (fixed, r.data[1:XAM2.aux_offset(r)], collect(aux_data(r)::XAM2.XAMAuxData.BAM.Auxiliary))
end

cigar_string(r::XAMRecord) = String(MemoryView(CIGAR(XAM2.cigar(r))))

# Convert the SAM line to a BAM record, returning the record or the error
function try_sam_to_bam(line::String, refs::References = References(["chr1", "chr2"], [1000, 500]))
    return XAM2.try_copy!(BAMRecord(refs), try_parse_line(line)::SAMRecord)
end

@testset "SAM <-> BAM conversion" begin
    bam_records = read_bam(EXAMPLE_BAM)
    sam_records = SAMReader(collect, open(EXAMPLE_SAM))
    refs = first(bam_records).refs

    @testset "SAM to BAM reproduces samtools' records" begin
        record = BAMRecord(refs)
        for (sam, bam) in zip(sam_records, bam_records)
            @test bam_fields(copy!(record, sam)) == bam_fields(bam)
        end
        # Missing quality is all 0xff
        r2 = copy!(record, sam_records[3])
        @test quality_scores(r2) === nothing
        @test all(==(0xff), r2.data[(XAM2.qual_offset(r2) + 1):XAM2.aux_offset(r2)])
    end

    @testset "Errors are returned, and reset the record" begin
        @test error_kind(try_sam_to_bam(with_field(3, "chr3"))) == XAMErrors.unknown_reference
        @test error_kind(try_sam_to_bam(with_field(7, "chr3"))) == XAMErrors.unknown_reference
        @test error_kind(try_sam_to_bam(with_field(6, "1M1H2M"))) == XAMErrors.invalid_cigar
        @test error_kind(try_sam_to_bam(with_field(6, "5M"))) == XAMErrors.cigar_query_length_mismatch
        @test error_kind(try_sam_to_bam(sam_line([VALID_FIELDS; "AB:i:x"]))) == XAMErrors.invalid_aux_data
        @test error_kind(try_sam_to_bam(sam_line([VALID_FIELDS; "AB"]))) == XAMErrors.invalid_aux_data

        record = copy(bam_records[1])
        sam = try_parse_line(with_field(3, "chr9"))
        @test XAM2.try_copy!(record, sam) == XAMError(XAMErrors.unknown_reference)
        @test is_reset(record)
        @test_throws XAMError(XAMErrors.unknown_reference) copy!(record, sam)
        @test is_reset(record)
    end

    @testset "Sequence characters not in BAM become N" begin
        record = try_sam_to_bam(with_field(10, "A.zt"))
        @test string(sequence(record)) == "ANNT"
    end

    @testset "Position of first =" begin
        for seq in ["A", "=", "AC", "A=", "=A", "=A=", "AC=", "ACGTA=", "A"^300 * "=" * "A"^10, "A"^301]
            record = try_sam_to_bam(sam_line(["r", "4", "*", "0", "0", "*", "*", "0", "0", seq, "*"]))
            @test XAM2.first_eq_seq_pos(record) == findfirst(==('='), seq)
        end
        @test XAM2.first_eq_seq_pos(BAMRecord(References(String[], Int[]))) === nothing
    end

    @testset "Undefined flag bits are cleared" begin
        record = try_parse_line(with_field(2, "65535"))
        @test flags(record) == SAMFlags(0xffff)
        @test Integer(flags(record)) == 0x0fff
        @test length(flags(record)) == 12
    end

    @testset "Bin" begin
        for T in (UInt32, UInt64)
            @test XAM2.reg2bin(T(0), T(1)) === UInt16(4681)
            @test XAM2.reg2bin(T(0), T(1 << 14)) == 4681
            @test XAM2.reg2bin(T(1 << 14), T((1 << 14) + 1)) == 4682
            @test XAM2.reg2bin(T(0), T((1 << 14) + 1)) == 585
            @test XAM2.reg2bin(T(0), T(1 << 26)) == 1
            @test XAM2.reg2bin(T(0), T((1 << 26) + 1)) == 0
            # Beyond BAI's range, the bin is truncated, like htslib
            @test XAM2.reg2bin(T(1 << 30), T((1 << 30) + 1)) == (4681 + (1 << 16)) % UInt16
        end
        # Unmapped record with no position
        record = try_sam_to_bam(sam_line(["r", "4", "*", "0", "0", "*", "*", "0", "0", "A", "*"]))
        @test record.bin == 4680
        # Unmapped record with a position
        record = try_sam_to_bam(sam_line(["r", "4", "chr1", "5", "0", "*", "*", "0", "0", "A", "*"]))
        @test record.bin == XAM2.reg2bin(UInt32(4), UInt32(5))
        # Zero reference length is binned as length 1
        record = try_sam_to_bam(sam_line(["r", "0", "chr1", "1", "0", "2S", "*", "0", "0", "AC", "*"]))
        @test record.bin == 4681
    end

    @testset "CIGARs with more than 65535 operations" begin
        n_ops = 70_000
        cigar_str = "1M1I"^(n_ops ÷ 2)
        seq = "A"^n_ops
        refs = References(["chr1"], [100_000])
        line = sam_line(["r", "0", "chr1", "1", "0", cigar_str, "*", "0", "0", seq, "*", "CG:Z:x", "NM:i:1"])
        record = try_sam_to_bam(line, refs)
        @test record isa BAMRecord
        # Placeholder 70000S35000N, and the CIGAR in the CG tag replacing the existing one
        @test record.n_cigar_op == 2
        @test cigar_string(record) == cigar_str
        @test mapping_span(record) == 1:(n_ops ÷ 2)
        @test count(t -> first(t) == AuxTag("CG"), aux_data(record)) == 1
        @test aux_data(record)["NM"] == 1

        # Survives writing and reading
        written = only(read_bam(write_bam(mktempdir(), [record], parse_header("@SQ\tSN:chr1\tLN:100000"), refs)))
        @test cigar_string(written) == cigar_str

        # Converting back puts the CIGAR in the CIGAR field, and removes the CG tag
        sam = copy!(SAMRecord(), written)
        @test cigar_string(sam) == cigar_str
        @test !haskey(aux_data(sam), "CG")
        @test aux_data(sam)["NM"] == 1
        @test String(copy(sam.data)) == replace(line, "\tCG:Z:x" => "")

        # Also when CG was the only tag: no trailing tab
        line = sam_line(["r", "0", "chr1", "1", "0", cigar_str, "*", "0", "0", seq, "*"])
        sam = copy!(SAMRecord(), try_sam_to_bam(line, refs))
        @test String(copy(sam.data)) == line
        @test isempty(aux_data(sam))
    end

    @testset "BAM to SAM errors" begin
        # Malformed auxiliary data
        record = copy(bam_records[1])
        push!(record.data, UInt8('X'))
        sam = copy!(SAMRecord(), bam_records[2])
        @test XAM2.try_copy!(sam, record) == XAMError(XAMErrors.invalid_aux_data)
        @test is_reset(sam)
    end
end
