using Test
using XAM2
using XAM2: XAMErrors, Reference, try_write, try_parse_header, seek_first_record
using BGZFLib: BGZFReader, BGZFWriter, BGZFError, VirtualOffset
using MemoryViews: ImmutableMemoryView, MemoryView
using CIGARStrings: CIGARStrings, CIGAR, CIGARError
using XAM2: AuxTag

const DATADIR = joinpath(@__DIR__, "data")
# Created from example.sam with `samtools view --no-PG -b -o example.bam example.sam`
const EXAMPLE_BAM = joinpath(DATADIR, "example.bam")
const EXAMPLE_SAM = joinpath(DATADIR, "example.sam")

parse_header(s::String) = try_parse_header(ImmutableMemoryView(s))::SAMHeader

# The kind of an error value, or the value itself if it is not an error
error_kind(x) = x isa XAMError ? x.kind : x

read_bam(path::String) = BAMReader(collect, open(path); n_workers = 1)

decompress(path::String) = open(io -> read(BGZFReader(io; n_workers = 1)), path)

# Write `records` to a new BAM file in `dir`, and return its path
function write_bam(dir::String, records, args...)
    path = tempname(dir)
    BAMWriter(open(path, "w"), args...; n_workers = 1) do writer
        for record in records
            write(writer, record)
        end
    end
    return path
end

summarize(r::BAMRecord) = (
    query_name(r),
    Integer(flags(r)),
    reference_name(r),
    left_position(r),
    next_reference_name(r),
    sequence(r) === nothing ? nothing : string(sequence(r)),
)

@testset "References" begin
    refs = References(["chr1", "chr2"], [1000, 500])
    @test length(refs) == 2
    @test refs[1] isa Reference
    @test refs[2].name == "chr2"
    @test refs[2].length == 500
    @test [(r.name, r.length) for r in refs] == [("chr1", 1000), ("chr2", 500)]
    @test refs[1] == References(["x", "chr1"], [1, 1000])[2]
    @test refs[1] != References(["chr1"], [1001])[1]
    @test refs == References(["chr1", "chr2"], [1000, 500])
    @test refs != References(["chr1", "chr2"], [1000, 501])
    @test sprint(show, refs[1]) == "Reference(\"chr1\", 1000)"
    @test_throws BoundsError refs[3]
    @test isempty(References(String[], Int[]))

    # The name index is built lazily, on the first lookup
    has_index(refs) = refs.lazy_index.index !== nothing
    @test !has_index(refs)
    @test XAM2.reference_id(refs, "chr2") == 1
    @test XAM2.reference_id(refs, "chr3") === nothing
    @test has_index(refs)

    # Invalid references
    @test_throws ArgumentError References(["a", "a"], [1, 2])
    @test_throws ArgumentError References([""], [1])
    @test_throws ArgumentError References(["a\0b"], [1])
    @test_throws ArgumentError References(["a"], [0])
    @test_throws ArgumentError References(["a"], [Int64(typemax(Int32)) + 1])
    @test_throws ArgumentError References(["a", "b"], [1])

    # From header
    header = parse_header("@HD\tVN:1.6\n@SQ\tSN:x\tLN:10\n@SQ\tSN:y\tLN:20\n")
    @test References(header) == References(["x", "y"], [10, 20])
    @test isempty(References(parse_header("@HD\tVN:1.6\n")))
    # Duplicate names are rejected by the header parser
    @test try_parse_header(ImmutableMemoryView("@SQ\tSN:x\tLN:10\n@SQ\tSN:x\tLN:20\n")) isa XAMError

    # Validating a header against references
    mismatch = XAMError(XAMErrors.header_reference_mismatch)
    @test XAM2.try_validate(header, References(["x", "y"], [10, 20])) === nothing
    @test XAM2.try_validate(parse_header("@HD\tVN:1.6\n"), References(["x"], [10])) === nothing
    @test XAM2.try_validate(header, References(["y", "x"], [20, 10])) == mismatch
    @test XAM2.try_validate(header, References(["x", "y"], [10, 21])) == mismatch
    @test XAM2.try_validate(header, References(["x"], [10])) == mismatch
end

@testset "BAMReader" begin
    reader = BAMReader(open(EXAMPLE_BAM); n_workers = 1)
    @test Base.IteratorSize(typeof(reader)) == Base.SizeUnknown()
    @test eltype(reader) == BAMRecord

    header = get_header(reader)
    @test length(header.seqs) == 3
    @test length(header.read_groups) == 1
    @test length(header.programs) == 1
    @test length(header.comments) == 1
    @test references(reader) == References(["chr1", "chr2", "chr3"], [1000, 500, 2000])

    records = collect(reader)
    @test records isa Vector{BAMRecord}
    @test map(summarize, records) == [
        ("r1", 0x0063, "chr1", 5, "chr1", "ACGT"),
        ("r1", 0x0093, "chr1", 20, "chr1", "TTTA"),
        ("r2", 0x0041, "chr1", 100, "chr2", "ACGTA"),
        ("r2", 0x0081, "chr2", 50, "chr1", "GGGCC"),
        ("r3", 0x0000, "chr3", 10, nothing, "ACGTA"),
        ("r4", 0x0004, nothing, nothing, nothing, "ACG"),
    ]
    @test map(is_mapped, records) == [true, true, true, true, true, false]
    @test map(mapping_span, records) == [5:8, 20:22, 100:102, 50:54, 10:15, nothing]
    @test mapping_span(records[2], XAM2.cigar(records[2])) == 20:22
    @test map(r -> String(MemoryView(CIGAR(XAM2.cigar(r)))), records[1:5]) == ["4M", "2M1I1M", "3M2S", "5M", "3M1D2M"]
    @test XAM2.cigar(records[6]) === nothing
    @test map(reference_length, records) == [1000, 1000, 1000, 500, 2000, nothing]
    @test map(reference_index, records) == [1, 1, 1, 2, 3, nothing]
    @test map(next_reference_index, records) == [1, 1, 2, 1, nothing, nothing]
    @test map(next_position, records) == [20, 5, 50, 100, nothing, nothing]
    @test map(mapping_quality, records) == [60, 60, 30, 30, nothing, 0]
    @test map(template_length, records) == [19, -19, nothing, nothing, nothing, nothing]
    @test all(r -> references(r) === references(reader), records)
    for r in records[1:5]
        @test references(r)[reference_index(r)] == XAM2.reference(r)
    end
    for r in records[1:4]
        @test references(r)[next_reference_index(r)].name == next_reference_name(r)
    end
    @test XAM2.reference(records[6]) === nothing
    @test quality_scores(records[2]) == [0x02, 0x03, 0x04, 0x05]
    @test quality_scores(records[3]) === nothing
    @test String(aux_data(records[5])["XA"]) == "hello"
    @test isempty(aux_data(records[2]))

    # Copying a record gives an independent copy
    original = records[2]
    @test XAM2.cigar(original) isa XAM2.BAMCIGAR
    cp = copy(original)
    @test XAM2.cigar(cp) == XAM2.cigar(original)
    # Each CIGAR aliases its own record's data
    @test pointer(MemoryView(XAM2.cigar(cp))) == pointer(cp.data) + cp.l_read_name
    @test pointer(MemoryView(XAM2.cigar(original))) == pointer(original.data) + original.l_read_name

    # Reading does not build the name index of the references
    @test references(reader).lazy_index.index === nothing

    # All records share the reader's references
    @test all(r -> r.refs === references(reader), records)

    # Seeking back reads the same records again
    seek_first_record(reader)
    @test map(summarize, collect(reader)) == map(summarize, records)
    close(reader)

    # With copy = false, the same record is returned every time
    BAMReader(open(EXAMPLE_BAM); n_workers = 1, copy = false) do reader
        record = XAM2.try_read!(reader)
        @test record isa BAMRecord
        @test XAM2.try_read!(reader) === record
    end
end

@testset "BAMWriter" begin
    reader = BAMReader(open(EXAMPLE_BAM); n_workers = 1)
    header = get_header(reader)
    refs = references(reader)
    records = collect(reader)
    close(reader)
    expected = map(summarize, records)
    no_sq_header = parse_header("@HD\tVN:1.6\tSO:unsorted\n")

    mktempdir() do dir
        @testset "Raw copy" begin
            # Shared references: Output is byte-identical to the input
            path = write_bam(dir, records, header, refs)
            @test decompress(path) == decompress(EXAMPLE_BAM)

            # References created from the header are content-identical
            path = write_bam(dir, records, header)
            @test decompress(path) == decompress(EXAMPLE_BAM)
            equal_refs = References(["chr1", "chr2", "chr3"], [1000, 500, 2000])
            path = write_bam(dir, records, header, equal_refs)
            @test decompress(path) == decompress(EXAMPLE_BAM)
            # No name lookups were needed
            @test equal_refs.lazy_index.index === nothing
        end

        @testset "Reference translation" begin
            # Reordered references
            new_refs = References(["chr3", "chr1", "chr2"], [2000, 1000, 500])
            path = write_bam(dir, records, no_sq_header, new_refs)
            @test new_refs.lazy_index.index !== nothing
            result = read_bam(path)
            @test map(summarize, result) == expected
            @test references(BAMReader(open(path); n_workers = 1)) == new_refs
            @test all(r -> r.refs !== refs, result)

            # Records from different sources, alternating, to the same writer
            refs_b = References(["chr2", "chr3", "chr1"], [500, 2000, 1000])
            records_b = read_bam(write_bam(dir, records, no_sq_header, refs_b))
            interleaved = collect(Iterators.flatten(zip(records, records_b)))
            path = write_bam(dir, interleaved, header, refs)
            @test map(summarize, read_bam(path)) == map(summarize, interleaved)
        end

        @testset "Missing or mismatched references" begin
            path = tempname(dir)
            BAMWriter(open(path, "w"), no_sq_header, References(["chr1"], [1000]); n_workers = 1) do writer
                @test map(r -> error_kind(try_write(writer, r)), records) == [
                    60, 57,
                    XAMErrors.unknown_reference, # mate on chr2
                    XAMErrors.unknown_reference, # chr2
                    XAMErrors.unknown_reference, # chr3
                    44,
                ]
                @test_throws XAMError write(writer, records[5])
            end
            # Failed writes write nothing
            @test map(summarize, read_bam(path)) == expected[[1, 2, 6]]

            BAMWriter(IOBuffer(), no_sq_header, References(["chr1"], [999]); n_workers = 1) do writer
                @test try_write(writer, records[1]) == XAMError(XAMErrors.reference_length_mismatch)
                @test try_write(writer, records[6]) isa Int
            end
        end

        @testset "Empty record" begin
            path = tempname(dir)
            BAMWriter(open(path, "w"), header, refs; n_workers = 1) do writer
                write(writer, BAMRecord(writer))
            end
            (record,) = read_bam(path)
            @test !is_mapped(record)
            @test query_name(record) === nothing
            @test reference_name(record) === nothing
        end

        @testset "Constructor and close" begin
            # @SQ lines in header must match references
            @test_throws ArgumentError BAMWriter(IOBuffer(), header, References(["chr1"], [1000]))
            @test_throws ArgumentError BAMWriter(
                IOBuffer(), header, References(["chr1", "chr2", "chr3"], [1000, 500, 2001])
            )

            writer = BAMWriter(IOBuffer(), header; n_workers = 1)
            @test references(writer) == refs
            close(writer)
            close(writer) # idempotent
            @test_throws Exception write(writer, records[1])
        end
    end
end

include("malformed.jl")
include("header.jl")
include("sam.jl")
include("convert.jl")
include("qualityview.jl")
include("biosequences.jl")
