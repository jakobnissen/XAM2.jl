# Tests of the BioSequences extension

using BioSequences: LongDNA, @dna_str

@testset "BioSequences extension" begin
    refs = References(["chr1"], [10])
    reference = dna"ACGTACGTAC"
    function record(rname, pos, cigar, seq)
        line = sam_line(["r", "0", rname, pos, "60", cigar, "*", "0", "0", seq, "*"])
        return try_sam_to_bam(line, refs)::BAMRecord
    end
    mapped(pos, cigar, seq) = record("chr1", pos, cigar, seq)
    seq(r) = XAM2.sequence(LongDNA{4}, r, reference)

    @testset "try_sequence" begin
        @test XAM2.try_sequence(LongDNA{4}, mapped("1", "*", "*")) === nothing
        @test XAM2.try_sequence(LongDNA{4}, mapped("1", "*", "A=")) === XAM2.HasEq()
        # Lengths around the chunk size of 16 bases, and odd lengths with a padding nibble
        for len in [1, 2, 15, 16, 17, 33]
            s = ("ACGTMRSVWYHKDBN"^3)[1:len]
            dna = XAM2.try_sequence(LongDNA{4}, mapped("1", "*", s))
            @test dna == LongDNA{4}(s)
            @test hash(dna) == hash(LongDNA{4}(s))
        end
    end

    @testset "sequence with reference" begin
        @test seq(mapped("1", "*", "*")) === nothing
        @test seq(mapped("1", "4M", "ACGA")) == dna"ACGA"
        @test seq(mapped("3", "4M", "=T=A")) == dna"GTAA"
        # Soft clip, then a deletion
        @test seq(mapped("2", "1S2M1D2M", "A=G=C")) == dna"ACGAC"
        @test seq(mapped("7", "4M", "====")) == dna"GTAC"

        # Missing reference, position or CIGAR
        @test seq(record("*", "0", "*", "=A")) === nothing
        @test seq(mapped("0", "*", "=A")) === nothing
        @test seq(mapped("1", "*", "=A")) === nothing

        # `=` not aligned to a reference base
        @test error_kind(seq(mapped("1", "1S1M", "=A"))) == XAMErrors.invalid_sequence
        @test error_kind(seq(mapped("1", "1M1I1M", "A=A"))) == XAMErrors.invalid_sequence
        # Alignment extends past the end of the reference
        @test error_kind(seq(mapped("8", "4M", "===="))) == XAMErrors.invalid_position
        @test error_kind(seq(mapped("8", "2M2D", "=="))) == XAMErrors.invalid_position

        # Malformed CIGAR, with an invalid operation
        r = mapped("1", "2M", "=A")
        r.data[XAM2.cigar_offset(r) + 1] |= 0x0f
        @test error_kind(XAM2.cigar(r)) == XAMErrors.invalid_cigar
        @test XAM2.cigar(r).wrapped_error isa CIGARError
        @test error_kind(seq(r)) == XAMErrors.invalid_cigar

        @test_throws ArgumentError XAM2.sequence(LongDNA{4}, mapped("1", "2M", "=A"), dna"ACGT")
    end
end
