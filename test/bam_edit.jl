@testset "Fixed-field setters" begin
    refs = References(["chr1", "chr2"], [1000, 500])
    header = parse_header("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500")
    record = BAMRecord(refs)

    # Setters return the record, and getters return what was set
    flags = SAMFlags([SAMFlag.PAIRED, SAMFlag.READ1])
    @test set_flags!(record, flags) === record
    @test XAM2.flags(record) == flags

    for (set!, get, values) in [
            (set_mapping_quality!, mapping_quality, [0, 60, 254, nothing]),
            (set_left_position!, left_position, [1, 17, typemax(Int32), nothing]),
            (set_next_position!, next_position, [1, 17, typemax(Int32), nothing]),
            (set_reference!, reference_index, [1, 2, nothing]),
            (set_next_reference!, next_reference_index, [1, 2, nothing]),
            (set_template_length!, template_length, [1, -1, typemin(Int32), typemax(Int32), nothing]),
        ]
        for value in values
            @test set!(record, value) === record
            @test get(record) == value
        end
    end

    # References by name
    set_reference!(record, "chr2")
    @test reference_name(record) == "chr2"
    set_next_reference!(record, "chr1")
    @test next_reference_name(record) == "chr1"

    # Invalid input throws, and leaves the record unchanged
    before = copy(record)
    @test_throws ArgumentError set_mapping_quality!(record, 255)
    @test_throws ArgumentError set_mapping_quality!(record, -1)
    @test_throws ArgumentError set_left_position!(record, 0)
    @test_throws ArgumentError set_left_position!(record, typemax(Int32) + 1)
    @test_throws ArgumentError set_next_position!(record, -1)
    @test_throws BoundsError set_reference!(record, 0)
    @test_throws BoundsError set_reference!(record, 3)
    @test_throws KeyError set_reference!(record, "chr3")
    @test_throws KeyError set_next_reference!(record, "chr3")
    @test_throws ArgumentError set_template_length!(record, 0)
    @test_throws ArgumentError set_template_length!(record, typemax(Int32) + 1)
    @test record == before

    # An edited record survives writing and reading
    set_flags!(record, SAMFlags([SAMFlag.PAIRED]))
    set_reference!(record, "chr1")
    set_left_position!(record, 100)
    set_mapping_quality!(record, 30)
    set_next_reference!(record, 2)
    set_next_position!(record, 50)
    set_template_length!(record, -250)
    written = only(read_bam(write_bam(mktempdir(), [record], header, refs)))
    @test written == record
end
