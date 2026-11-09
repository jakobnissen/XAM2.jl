# Everything in a BAM record but the name, as comparable values that don't alias the record
function content(r::BAMRecord)
    cig = XAM2.cigar(r)
    seq = sequence(r)
    qual = quality_scores(r)
    own(v) = v isa AbstractString ? String(v) : v isa AbstractVector ? collect(v) : v
    return (
        XAM2.fixed_fields(r)[[1, 2, 4:end...]],
        cig === nothing ? nothing : String(MemoryView(CIGAR(cig))),
        seq === nothing ? nothing : String(seq),
        qual === nothing ? nothing : collect(qual),
        [k => own(v) for (k, v) in aux_data(r)::XAM2.XAMAuxData.BAM.Auxiliary],
    )
end

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

@testset "Query name setter" begin
    # Everything but the name, as comparable values
    records = read_bam(EXAMPLE_BAM)
    header = BAMReader(get_header, open(EXAMPLE_BAM))
    names = ["x", "r1", "!?A~", "a"^254, nothing]
    for original in records, name in names
        record = copy(original)
        @test set_query_name!(record, name) === record
        @test query_name(record) == name
        @test content(record) == content(original)
        written = only(read_bam(write_bam(mktempdir(), [record], header, original.refs)))
        @test written == record
        # Setting the original name back restores the record
        set_query_name!(record, query_name(original) === nothing ? nothing : String(query_name(original)))
        @test record == original
    end

    # Invalid names throw, and leave the record unchanged
    record = copy(first(records))
    for name in ["", "*", "a@b", "a b", "a\0b", "a\tb", "é", "a"^255]
        @test_throws ArgumentError set_query_name!(record, name)
        @test record == first(records)
    end
    # "*" is only special on its own
    @test query_name(set_query_name!(record, "**")) == "**"
end

@testset "Query name setter with aliased input" begin
    refs = References(["chr1"], [1000])
    function make_record(name, umi)
        fields = [name, "0", "chr1", "5", "60", "4M", "*", "0", "0", "ACGT", "IIII", "RX:Z:$umi", "NM:i:1"]
        return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
    end

    # The record's own name, and substrings of it, which shrink the name
    for (name, new) in [
            ("abcdef", r -> query_name(r)),
            ("abcdef", r -> SubString(query_name(r), 4)),
            ("abcdef", r -> SubString(query_name(r), 2, 3)),
            ("abcdef", r -> SubString(query_name(r), 1, 1)),
        ]
        record = make_record(name, "UMI")
        expected = String(new(record))
        before = content(record)
        set_query_name!(record, new(record))
        @test query_name(record) == expected
        @test content(record) == before
    end

    # Auxiliary string values, which are after the name, and grow or shrink it
    for (name, umi) in [("r", "ACGTACGTACGTACGT"), ("a_long_read_name", "ACG"), ("rrr", "AAA")]
        for new in [r -> aux_data(r)["RX"], r -> SubString(aux_data(r)["RX"], 2)]
            record = make_record(name, umi)
            expected = String(new(record))
            before = content(record)
            set_query_name!(record, new(record))
            @test query_name(record) == expected
            @test content(record) == before
            @test aux_data(record)["RX"] == umi
        end
    end

    # Only views of the record's own data are detected
    record = make_record("abc", "UMI")
    other = make_record("abc", "UMI")
    data = record.data
    indices(x) = XAM2.data_indices(data, x)
    @test indices(codeunits(query_name(record))) == 1:3
    @test indices(codeunits(SubString(query_name(record), 2))) == 2:3
    @test indices(codeunits(query_name(other))) === nothing
    @test indices(codeunits("abc")) === nothing
    @test indices(codeunits(SubString(query_name(record), 2, 1))) === nothing
    # Contiguous SubArrays
    @test indices(view(data, 2:3)) == 2:3
    @test indices(view(data, :)) == eachindex(data)
    @test indices(SubArray(MemoryView(data), (2:3,))) == 2:3
    @test indices(view(other.data, 2:3)) === nothing
    # Non-contiguous SubArrays are not checked
    @test indices(view(data, 1:2:5)) === nothing
    @test indices(view(data, [1, 2])) === nothing
end

# A string with UInt16 code units, like UTF-16, supporting only the BMP
struct UInt16String <: AbstractString
    units::Vector{UInt16}
end
UInt16String(s::String) = UInt16String([UInt16(c) for c in s])
Base.ncodeunits(s::UInt16String) = length(s.units)
Base.codeunit(::UInt16String) = UInt16
Base.codeunit(s::UInt16String, i::Integer) = s.units[i]
Base.isvalid(s::UInt16String, i::Int) = checkbounds(Bool, s.units, i)
Base.iterate(s::UInt16String, i::Int = 1) = i > length(s.units) ? nothing : (Char(s.units[i]), i + 1)

# A string whose code units are all in the ASCII range, but which encodes a different
# string, like UTF-7, which encodes "a&b" as "a+ACY-b".
struct EncodedString <: AbstractString
    units::Vector{UInt8}
    chars::Vector{Char}
end
Base.ncodeunits(s::EncodedString) = length(s.units)
Base.codeunit(::EncodedString) = UInt8
Base.codeunit(s::EncodedString, i::Integer) = s.units[i]
Base.isvalid(s::EncodedString, i::Int) = i == 1 || i == length(s.units) + 1
function Base.iterate(s::EncodedString, i::Int = 1)
    return i > length(s.chars) ? nothing : (s.chars[i], i + 1)
end

@testset "Query name setter with non-UTF-8 strings" begin
    record = BAMRecord(References(String[], Int[]))
    @test query_name(set_query_name!(record, UInt16String("read_1"))) == "read_1"
    @test_throws ArgumentError set_query_name!(record, UInt16String("a b"))
    @test query_name(record) == "read_1"
    # The name is set from the characters, not the code units
    @test query_name(set_query_name!(record, EncodedString(collect(b"a+ACY-b"), collect("a&b")))) == "a&b"
    @test query_name(set_query_name!(record, Test.GenericString("read_2"))) == "read_2"
end
