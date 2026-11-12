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
    flags = SAMFlags([SAMFlag.PAIRED, SAMFlag.READ1, SAMFlag.UNMAP])
    @test set_flags!(record, flags) === record
    @test XAM2.flags(record) == flags

    for (set!, get, values) in [
            (set_mapping_quality!, mapping_quality, [0, 60, 254, nothing]),
            (set_template_length!, template_length, [1, -1, typemin(Int32), typemax(Int32), nothing]),
        ]
        for value in values
            @test set!(record, value) === record
            @test get(record) == value
        end
    end

    # Placement, by index or name
    for (set!, get_ref, get_pos) in [
            (set_placement!, reference_index, left_position),
            (set_next_placement!, next_reference_index, next_position),
        ]
        for (ref, pos) in [(1, 1), (2, 17), ("chr1", 1000), ("chr2", 500)]
            @test set!(record, ref, pos) === record
            @test get_ref(record) == (ref isa String ? (ref == "chr1" ? 1 : 2) : ref)
            @test get_pos(record) == pos
        end
        @test set!(record, nothing) === record
        @test get_ref(record) === nothing
        @test get_pos(record) === nothing
    end

    # Invalid input throws, and leaves the record unchanged
    set_placement!(record, 1, 10)
    before = copy(record)
    @test_throws ArgumentError set_mapping_quality!(record, 255)
    @test_throws ArgumentError set_mapping_quality!(record, -1)
    for set! in [set_placement!, set_next_placement!]
        @test_throws ArgumentError set!(record, 1, 0)
        @test_throws ArgumentError set!(record, 1, typemax(Int32) + 1)
        @test_throws ArgumentError set!(record, 1, 1001)
        @test_throws ArgumentError set!(record, "chr2", 501)
        @test_throws BoundsError set!(record, 0, 1)
        @test_throws BoundsError set!(record, 3, 1)
        @test_throws KeyError set!(record, "chr3", 1)
    end
    @test_throws ArgumentError set_template_length!(record, 0)
    @test_throws ArgumentError set_template_length!(record, typemax(Int32) + 1)
    @test record == before

    # A mapped record can't be unplaced, and an unplaced record can't be mapped
    mapped = SAMFlags([SAMFlag.PAIRED])
    set_flags!(record, mapped)
    @test_throws ArgumentError set_placement!(record, nothing)
    set_flags!(record, flags)
    set_placement!(record, nothing)
    before = copy(record)
    @test_throws ArgumentError set_flags!(record, mapped)
    @test record == before
    # Other flags can be changed on an unplaced record
    set_flags!(record, SAMFlags([SAMFlag.UNMAP]))

    # A record that is already mapped and unplaced can't keep that state
    broken = try_sam_to_bam(sam_line(["r", "0", "*", "0", "255", "*", "*", "0", "0", "*", "*"]), refs)::BAMRecord
    @test_throws ArgumentError set_flags!(broken, SAMFlags([SAMFlag.DUP]))
    @test set_flags!(broken, SAMFlags([SAMFlag.UNMAP])) === broken

    # An edited record survives writing and reading
    set_placement!(record, "chr1", 100)
    set_flags!(record, mapped)
    set_mapping_quality!(record, 30)
    set_next_placement!(record, 2, 50)
    set_template_length!(record, -250)
    written = only(read_bam(write_bam(mktempdir(), [record], header, refs)))
    @test written == record
end

@testset "CIGAR setter" begin
    refs = References(["chr1"], [10^9])
    header = parse_header("@SQ\tSN:chr1\tLN:1000000000")
    function make_record(cigar, seq; aux = ["NM:i:1"])
        fields = ["r", "0", "chr1", "5", "60", cigar, "*", "0", "0", seq, "*", aux...]
        return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
    end
    roundtrip(r) = only(read_bam(write_bam(mktempdir(), [r], header, refs)))
    bam(s) = BAMCIGAR(CIGAR(s))

    # Grow, shrink, keep length, remove, and add
    record = make_record("4M", "ACGT")
    for new in ["2M1I1M", "1S2M2D1M", "4M", nothing, "4=", "1X3M"]
        cig = new === nothing ? nothing : bam(new)
        @test set_cigar!(record, cig) === record
        @test XAM2.cigar(record) == cig
        @test sequence(record) == "ACGT"
        @test aux_data(record)["NM"] == 1
        @test roundtrip(record) == record
    end

    # A record without a sequence takes any CIGAR
    record = make_record("*", "*")
    set_cigar!(record, bam("10M"))
    @test XAM2.cigar(record) == bam("10M")

    # Invalid input throws, and leaves the record unchanged
    record = make_record("4M", "ACGT")
    before = copy(record)
    @test_throws ArgumentError set_cigar!(record, bam("5M"))
    @test_throws ArgumentError set_cigar!(record, BAMCIGAR(UInt8[]))
    @test record == before

    # Aliased input: the record's own CIGAR, and a CIGAR stored in an aux array
    record = make_record("2M1I1M", "ACGT")
    set_cigar!(record, XAM2.cigar(record))
    @test XAM2.cigar(record) == bam("2M1I1M")
    @test record == make_record("2M1I1M", "ACGT")
    xc = join(reinterpret(UInt32, collect(MemoryView(bam("1M1I1M1I")))), ',')
    record = make_record("4M", "ACGT"; aux = ["XC:B:I,$xc"])
    set_cigar!(record, BAMCIGAR(parent(aux_data(record)["XC"])))
    @test XAM2.cigar(record) == bam("1M1I1M1I")
    @test roundtrip(record) == record

    # A stale CG tag is removed
    record = make_record("4M", "ACGT"; aux = ["CG:B:I,64", "NM:i:1"])
    set_cigar!(record, bam("4M"))
    @test !haskey(aux_data(record), "CG")
    @test aux_data(record)["NM"] == 1

    # Long CIGARs are stored in the CG tag
    n = 70_000
    long = bam("1M1I"^n)
    record = make_record("*", "A"^(2n))
    set_cigar!(record, long)
    @test XAM2.cigar(record) == long
    @test record.n_cigar_op == 2
    @test aux_data(record)["NM"] == 1
    @test roundtrip(record) == record
    # Setting the CIGAR from the CG tag itself, and replacing it with another long CIGAR
    set_cigar!(record, XAM2.cigar(record))
    @test XAM2.cigar(record) == long
    long2 = bam("2M2I"^(n ÷ 2))
    set_cigar!(record, long2)
    @test XAM2.cigar(record) == long2
    @test length(collect(aux_data(record))) == 2
    # The placeholder makes unplacing throw
    set_flags!(record, SAMFlags([SAMFlag.UNMAP]))
    @test_throws ArgumentError set_placement!(record, nothing)
    # Replacing with a short CIGAR removes the CG tag
    set_cigar!(record, bam("$(2n)M"))
    @test XAM2.cigar(record) == bam("$(2n)M")
    @test !haskey(aux_data(record), "CG")
    set_placement!(record, nothing)
    # Long CIGARs need a reference
    @test_throws ArgumentError set_cigar!(record, long)
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

@testset "Sequence setters" begin
    refs = References(["chr1"], [10^9])
    header = parse_header("@SQ\tSN:chr1\tLN:1000000000")
    function make_record(cigar, seq, qual = "*"; aux = ["NM:i:1"])
        fields = ["r", "0", "chr1", "5", "60", cigar, "*", "0", "0", seq, qual, aux...]
        return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
    end
    roundtrip(r) = only(read_bam(write_bam(mktempdir(), [r], header, refs)))
    bam(s) = BAMCIGAR(CIGAR(s))
    phred(s) = s == "*" ? nothing : [UInt8(c) - 0x21 for c in s]

    # Without a CIGAR, any length can be set. Compare with records built from SAM.
    cases = [("A", "I"), ("AC", "*"), ("ACG", "!#I"), ("ACGTN", "*"), ("*", "*"), ("=ACGTRYKMSWBDHVN", "!"^16)]
    for (seq0, qual0) in cases, (seq, qual) in cases
        record = make_record("*", seq0, qual0)
        s = seq == "*" ? nothing : seq
        @test set_sequence!(record, s, phred(qual)) === record
        @test record == make_record("*", seq, qual)
        @test roundtrip(record) == record
    end

    # Lowercase is stored as uppercase
    record = make_record("*", "*")
    set_sequence!(record, "acgT=mrsvwyhkdbn", nothing)
    @test sequence(record) == "ACGT=MRSVWYHKDBN"

    # A sequence from another record, also with odd length, and non-String strings
    for seq in ["ACGTA", "ACGT"]
        record = make_record("*", "*")
        set_sequence!(record, sequence(make_record("*", seq)), nothing)
        @test sequence(record) == seq
        set_sequence!(record, SubString(seq, 2), nothing)
        @test sequence(record) == seq[2:end]
        set_sequence!(record, SubString(sequence(make_record("*", seq)), 2), nothing)
        @test sequence(record) == seq[2:end]
    end

    # Invalid input throws, and leaves the record unchanged
    record = make_record("3M", "ACG", "III")
    before = copy(record)
    for (seq, qual) in [
            ("", nothing), ("AC GT", nothing), ("ACé", nothing), ("AC1", nothing),
            # Letters that are not IUPAC codes, and `.`, which are accepted when parsing SAM
            ("ACX", nothing), ("ACx", nothing), ("AC.", nothing), ("ACU", nothing),
            (nothing, UInt8[]), ("ACG", UInt8[1, 2]), ("ACG", UInt8[1, 2, 94]),
            # The length of a record with a CIGAR can't change
            ("ACGT", nothing), (nothing, nothing),
        ]
        @test_throws ArgumentError set_sequence!(record, seq, qual)
        @test record == before
    end
    # The same length is fine
    set_sequence!(record, "TTT", UInt8[93, 0, 1])
    @test record == make_record("3M", "TTT", "~!\"")
    # A record with a CIGAR, but no sequence, can't get a sequence
    record = make_record("3M", "*")
    @test_throws ArgumentError set_sequence!(record, "ACG", nothing)

    # Aliased input: the record's own sequence and quality, and an aux array
    record = make_record("*", "ACGTA", "ABCDE")
    set_sequence!(record, sequence(record), quality_scores(record))
    @test record == make_record("*", "ACGTA", "ABCDE")
    set_sequence!(record, SubString(sequence(record), 2), quality_scores(record)[2:end])
    @test record == make_record("*", "CGTA", "BCDE")
    record = make_record("*", "A", "I"; aux = ["XQ:B:C,1,2,3,4,5,6,7", "NM:i:1"])
    set_sequence!(record, "ACGTACG", aux_data(record)["XQ"])
    @test record == make_record("*", "ACGTACG", "\"#\$%&'("; aux = ["XQ:B:C,1,2,3,4,5,6,7", "NM:i:1"])
    set_sequence!(record, "AC", aux_data(record)["XQ"][3:4])
    @test record == make_record("*", "AC", "\$%"; aux = ["XQ:B:C,1,2,3,4,5,6,7", "NM:i:1"])

    # Setting the CIGAR and sequence together, e.g. trimming
    record = make_record("2S3M1I", "ACGTAC", "ABCDEF")
    @test set_cigar_and_sequence!(record, bam("3M"), SubString(sequence(record), 3, 5), quality_scores(record)[3:5]) === record
    @test record == make_record("3M", "GTA", "CDE")
    @test roundtrip(record) == record
    set_cigar_and_sequence!(record, bam("1M2I5M"), "ACGTACGT", nothing)
    @test record == make_record("1M2I5M", "ACGTACGT")
    set_cigar_and_sequence!(record, bam("4M"), nothing, nothing)
    @test record == make_record("4M", "*")
    set_cigar_and_sequence!(record, nothing, nothing, nothing)
    @test record == make_record("*", "*")
    set_cigar_and_sequence!(record, nothing, "AC", nothing)
    @test record == make_record("*", "AC")

    # Invalid input throws, and leaves the record unchanged
    record = make_record("3M", "ACG", "III")
    before = copy(record)
    @test_throws ArgumentError set_cigar_and_sequence!(record, bam("4M"), "ACG", nothing)
    @test_throws ArgumentError set_cigar_and_sequence!(record, BAMCIGAR(UInt8[]), "ACG", nothing)
    @test_throws ArgumentError set_cigar_and_sequence!(record, bam("2M"), "AC", UInt8[1])
    @test record == before

    # Long CIGARs, and removal of the CG tag
    n = 70_000
    long = bam("1M1I"^n)
    record = make_record("3M", "ACG")
    set_cigar_and_sequence!(record, long, "A"^(2n), nothing)
    @test XAM2.cigar(record) == long
    @test sequence(record) == "A"^(2n)
    @test aux_data(record)["NM"] == 1
    @test roundtrip(record) == record
    set_cigar_and_sequence!(record, XAM2.cigar(record), "C"^(2n), nothing)
    @test XAM2.cigar(record) == long
    @test sequence(record) == "C"^(2n)
    set_cigar_and_sequence!(record, bam("3M"), "ACG", nothing)
    @test record == make_record("3M", "ACG")
end

@testset "In-place sequence and quality" begin
    refs = References(["chr1"], [1000])
    function make_record(seq, qual = "*")
        fields = ["r", "0", "chr1", "5", "60", "*", "*", "0", "0", seq, qual, "NM:i:1"]
        return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
    end

    # Passing the record's own scores to a setter is checked for aliasing
    record = make_record("ACGTA", "ABCDE")
    set_sequence!(record, "CGTA", quality_scores(record)[2:end])
    @test record == make_record("CGTA", "BCDE")
    record = make_record("ACGTA", "ABCDE")
    set_sequence!(record, "ACGTA", quality_scores(record))
    @test record == make_record("ACGTA", "ABCDE")
    record = make_record("ACGTA", "ABCDE")
    set_sequence!(record, "CGT", view(quality_scores(record), 2:4))
    @test record == make_record("CGT", "BCD")

    # Raw sequence bytes, with odd and even lengths
    @test XAM2.unsafe_sequence_bytes(make_record("*")) isa XAM2.MutableMemoryView{UInt8}
    @test isempty(XAM2.unsafe_sequence_bytes(make_record("*")))
    record = make_record("ACGTA")
    bytes = XAM2.unsafe_sequence_bytes(record)
    @test bytes isa XAM2.MutableMemoryView{UInt8}
    @test bytes == [0x12, 0x48, 0x10]
    bytes[2] = 0xf0
    @test sequence(record) == "ACN=A"
    @test sequence(record) isa BAMSeqString
    @test length(XAM2.unsafe_sequence_bytes(make_record("ACGT"))) == 2
end

@testset "Slicing the sequence" begin
    refs = References(["chr1"], [10^9])
    header = parse_header("@SQ\tSN:chr1\tLN:1000000000")
    function make_record(cigar, seq, qual = "*"; aux = ["NM:i:1"])
        fields = ["r", "0", "chr1", "5", "60", cigar, "*", "0", "0", seq, qual, aux...]
        return try_sam_to_bam(sam_line(fields), refs)::BAMRecord
    end
    roundtrip(r) = only(read_bam(write_bam(mktempdir(), [r], header, refs)))
    bam(s) = BAMCIGAR(CIGAR(s))

    # Exhaustively against a naive implementation, with and without quality,
    # and with a CIGAR that shrinks, keeps or grows the CIGAR field
    bases = "=ACMGRSVTWYHKDBN"
    for len in 1:9, has_qual in (false, true), idx in [i:j for i in 1:len for j in i:len]
        seq = String([bases[mod1(3i + len, 16)] for i in 1:len])
        qual = has_qual ? String([Char(0x21 + mod(5i, 94)) for i in 1:len]) : "*"
        n = length(idx)
        expected_qual = has_qual ? qual[idx] : "*"

        record = make_record("*", seq, qual)
        @test slice_sequence!(record, idx) === record
        @test record == make_record("*", seq[idx], expected_qual)

        for (old_cigar, new_cigar) in [
                ("$(len)M", "$(n)M"),
                ("$(len)=", nothing),
                ("$(len)M", "1H" * "1M1D"^(n - 1) * "1M1H"),
            ]
            record = make_record(old_cigar, seq, qual)
            cig = new_cigar === nothing ? nothing : bam(new_cigar)
            @test slice_sequence!(record, idx, cig) === record
            @test record == make_record(something(new_cigar, "*"), seq[idx], expected_qual)
        end
    end

    # Missing quality stays missing, and the result can be written
    record = make_record("2S5M1S", "ACGTACGT")
    slice_sequence!(record, 3:7, bam("2H5M1H"))
    @test record == make_record("2H5M1H", "GTACG")
    @test roundtrip(record) == record

    # Long CIGARs, and removal of the CG tag
    n = 70_000
    long = bam("1M1I"^n)
    record = make_record("3M", "ACG")
    set_cigar_and_sequence!(record, long, "AC"^n, nothing)
    slice_sequence!(record, 2:(2n), bam("1I" * "1M1I"^(n - 1)))
    @test sequence(record) == "C" * "AC"^(n - 1)
    @test aux_data(record)["NM"] == 1
    long2 = bam("1M1I"^(n - 1))
    slice_sequence!(record, 2:(2n - 1), long2)
    @test XAM2.cigar(record) == long2
    @test sequence(record) == "AC"^(n - 1)
    @test roundtrip(record) == record
    slice_sequence!(record, 1:3, bam("3M"))
    @test record == make_record("3M", "ACA")

    # Invalid input throws, and leaves the record unchanged
    record = make_record("3M", "ACG", "III")
    before = copy(record)
    @test_throws ArgumentError slice_sequence!(record, 1:2)
    @test_throws ArgumentError slice_sequence!(record, 2:1, nothing)
    @test_throws ArgumentError slice_sequence!(record, 0:2, nothing)
    @test_throws ArgumentError slice_sequence!(record, 2:4, nothing)
    @test_throws ArgumentError slice_sequence!(record, 1:2, bam("3M"))
    @test record == before
    # The whole sequence can be kept without a CIGAR argument
    slice_sequence!(record, 1:3)
    @test record == before
    @test_throws ArgumentError slice_sequence!(make_record("*", "*"), 1:0)
end
