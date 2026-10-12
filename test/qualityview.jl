# Tests of QualityView, and the quality accessors that return it

using XAM2: QualityView
using XAM2: StringView

@testset "Quality scores" begin
    @testset "QualityView" begin
        v = UInt8[0x00, 0x21, 0xff]
        q = QualityView(ImmutableMemoryView(v), 0x21)
        @test q isa AbstractVector{UInt8}
        @test size(q) == (3,)
        # Arithmetic wraps around
        @test q == [0x21, 0x42, 0x20]
        @test_throws BoundsError q[0]
        @test_throws BoundsError q[4]
        @test collect(q) isa Vector{UInt8}
        @test collect(q) == [0x21, 0x42, 0x20]

        # Copying to arrays of other types
        @test copyto!(zeros(UInt8, 3), q) == q
        @test copyto!(Memory{UInt8}(undef, 3), q) == q
        @test copyto!(MemoryView(zeros(UInt8, 3)), q) == q
        @test copyto!(zeros(UInt8, 5), 2, q, 2, 2) == [0x00, 0x42, 0x20, 0x00, 0x00]

        # copyto! checks bounds
        @test_throws BoundsError copyto!(zeros(UInt8, 2), q)
        @test_throws BoundsError copyto!(zeros(UInt8, 5), 4, q, 1, 3)
        @test_throws BoundsError copyto!(zeros(UInt8, 5), 1, q, 2, 3)
        @test_throws ArgumentError copyto!(zeros(UInt8, 5), 1, q, 1, -1)
        @test copyto!(UInt8[7], 5, q, 5, 0) == [7]
        # Destinations with other element types, and Cartesian indexing
        @test copyto!(zeros(Int, 3), q) == [0x21, 0x42, 0x20]
        @test copyto!(zeros(UInt8, 3, 2), 2, q, 1, 3) == [0x00 0x20; 0x21 0x00; 0x42 0x00]
        @test copyto!(PermutedDimsArray(zeros(UInt8, 3), (1,)), q) == [0x21, 0x42, 0x20]
        @test isempty(Test.detect_ambiguities(XAM2))
    end

    bam_records = read_bam(EXAMPLE_BAM)
    sam_records = SAMReader(collect, open(EXAMPLE_SAM))

    @testset "SAM" begin
        for (sam, bam) in zip(sam_records, bam_records)
            scores = quality_scores(sam)
            expected = quality_scores(bam)
            if expected === nothing
                @test scores === nothing
            else
                @test scores isa QualityView
                @test scores == expected
            end
        end
        @test quality_string(sam_records[1]) isa StringView{ImmutableMemoryView{UInt8}}
        @test quality_scores(SAMRecord()) === nothing
        @test collect(quality_scores(try_parse_line(with_field(11, "!+I~")))) == [0, 10, 40, 93]
    end

    @testset "BAM" begin
        for (sam, bam) in zip(sam_records, bam_records)
            str = quality_string(bam)
            expected = quality_string(sam)
            if expected === nothing
                @test str === nothing
            else
                @test str isa StringView{QualityView}
                @test str == expected
                @test String(str) == expected
            end
        end
        @test quality_string(BAMRecord(first(bam_records).refs)) === nothing
    end
end
