# Tests of QualityView, and the quality accessors and setters

using XAM2: QualityView, unsafe_quality_bytes
using MemoryViews: MutableMemoryView
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
        # Immutable
        @test_throws Base.CanonicalIndexError q[1] = 0x00

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
        # Missing sequence, or only missing quality
        @test quality_scores(SAMRecord()) === nothing
        @test quality_string(SAMRecord()) === nothing
        @test quality_scores(try_parse_line(with_field(11, "*"))) === nothing
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
        @test quality_scores(first(bam_records)) isa ImmutableMemoryView{UInt8}
        @test quality_string(BAMRecord(first(bam_records).refs)) === nothing
        @test quality_scores(BAMRecord(first(bam_records).refs)) === nothing
    end

    @testset "Mutation" begin
        refs = References(["chr1"], [1000])
        function make_record(seq, qual = "*")
            fields = ["r", "0", "chr1", "5", "60", "*", "*", "0", "0", seq, qual]
            return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
        end

        # Raw bytes, also when the quality is missing
        bytes = unsafe_quality_bytes(make_record("*"))
        @test bytes isa MutableMemoryView{UInt8}
        @test isempty(bytes)
        record = make_record("ACG")
        bytes = unsafe_quality_bytes(record)
        @test bytes isa MutableMemoryView{UInt8}
        @test bytes == [0xff, 0xff, 0xff]
        copyto!(bytes, [0x00, 0x0a, 0x28])
        @test record == make_record("ACG", "!+I")
        bytes[1] = 0xff
        @test quality_scores(record) === nothing

        # Checked setter
        record = make_record("ACG")
        @test set_quality_scores!(record, [0x00, 0x0a, 0x5d]) === record
        @test record == make_record("ACG", "!+~")
        set_quality_scores!(record, nothing)
        @test record == make_record("ACG")
        @test_throws ArgumentError set_quality_scores!(make_record("*"), [0x00])
        @test_throws ArgumentError set_quality_scores!(make_record("*"), nothing)
        record = make_record("ACG", "!+I")
        @test_throws ArgumentError set_quality_scores!(record, [0x00, 0x00])
        @test_throws ArgumentError set_quality_scores!(record, [0x00, 0x00, 0x5e])
        @test_throws ArgumentError set_quality_scores!(record, [0xff, 0x00, 0x00])
        @test record == make_record("ACG", "!+I")

        # Sources that alias the record, or have other types
        set_quality_scores!(record, quality_scores(record))
        @test record == make_record("ACG", "!+I")
        set_quality_scores!(record, quality_scores(make_record("ACG", "?+I")))
        @test record == make_record("ACG", "?+I")
        set_quality_scores!(record, quality_scores(try_parse_line(with_field(11, "!+I~")))[1:3])
        @test record == make_record("ACG", "!+I")
    end
end
