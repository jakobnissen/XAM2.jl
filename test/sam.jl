# Tests of parsing SAM records, and of SAMReader

const VALID_FIELDS = ["r1", "99", "chr1", "5", "60", "4M", "chr2", "20", "-19", "ACGT", "IIII"]

sam_line(fields) = join(fields, '\t')

# Parse `line` into a new record, returning the record or the error
try_parse_line(line::String) = XAM2.try_parse!(SAMRecord(), ImmutableMemoryView(line))

# The line with field `i` (1-based) replaced by `value`
with_field(i::Int, value::String) = sam_line([j == i ? value : f for (j, f) in enumerate(VALID_FIELDS)])

function is_reset(record::SAMRecord)
    return query_name(record) === nothing &&
        reference_name(record) === nothing &&
        next_reference_name(record) === nothing &&
        left_position(record) === nothing &&
        XAM2.cigar(record) === nothing &&
        sequence(record) === nothing &&
        quality_string(record) === nothing &&
        isempty(aux_data(record)) &&
        !is_mapped(record)
end

@testset "SAM records" begin
    @testset "Empty record" begin
        @test is_reset(SAMRecord())
        @test XAM2.try_parse!(SAMRecord(), ImmutableMemoryView(SAMRecord().data)) isa SAMRecord
    end

    @testset "Valid records" begin
        record = try_parse_line(sam_line(VALID_FIELDS) * "\tNM:i:0\tXA:Z:hello world")
        @test record isa SAMRecord
        @test query_name(record) == "r1"
        @test Integer(flags(record)) == 99
        @test reference_name(record) == "chr1"
        @test left_position(record) === 5
        @test mapping_quality(record) === 60
        @test String(MemoryView(XAM2.cigar(record))) == "4M"
        # Regression: RNEXT overwrote RNAME
        @test next_reference_name(record) == "chr2"
        @test next_position(record) === 20
        @test template_length(record) === -19
        @test sequence(record) == "ACGT"
        @test quality_string(record) == "IIII"
        @test collect(aux_data(record)) == [AuxTag("NM") => 0, AuxTag("XA") => "hello world"]
        @test mapping_span(record) == 5:8
        @test mapping_span(record, XAM2.cigar(record)) == 5:8
        @test collect(mapping_blocks(record)) == [5:8]

        # Spliced record with a deletion
        record = try_parse_line(with_field(6, "1S1M1D1M5N1M"))
        @test collect(mapping_blocks(record)) == [5:7, 13:13]
        @test collect(mapping_blocks(record, XAM2.cigar(record))) == [5:7, 13:13]
        @test collect(mapping_blocks(record, XAM2.cigar(record), CIGARStrings.query)) == [5:5, 7:7, 13:13]
        @test mapping_blocks(try_parse_line(with_field(4, "0"))) === nothing

        # Missing fields
        record = try_parse_line("*\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*")
        @test is_reset(record)
        @test left_position(record) === nothing
        @test next_position(record) === nothing
        @test mapping_quality(record) === 0
        @test template_length(record) === nothing

        # Special values
        record = try_parse_line(sam_line(["r", "0", "chr1", "1", "0", "*", "=", "0", "+5", "AC=.n", "*"]))
        @test next_reference_name(record) == "chr1"
        @test template_length(record) == 5
        @test sequence(record) == "AC=.n"
        @test quality_string(record) === nothing
        @test XAM2.cigar(record) === nothing
        @test isempty(aux_data(record))

        # Limits of fields
        record = try_parse_line(
            sam_line(
                [
                    "!"^254, "65535", "a*=", "2147483647", "255", "1M",
                    "*", "2147483647", "-2147483648", "A", "~",
                ]
            )
        )
        @test record isa SAMRecord
        @test length(query_name(record)) == 254
        @test reference_name(record) == "a*="
        @test left_position(record) == typemax(Int32)
        @test next_position(record) == typemax(Int32)
        @test mapping_quality(record) === nothing
        @test template_length(record) == typemin(Int32)
        @test template_length(try_parse_line(with_field(9, "2147483647"))) == typemax(Int32)

        # The CIGAR is validated lazily
        record = try_parse_line(with_field(6, "3M"))
        @test record isa SAMRecord
        @test XAM2.cigar(record) == XAMError(XAMErrors.cigar_query_length_mismatch)
        @test mapping_span(record) == XAMError(XAMErrors.cigar_query_length_mismatch)
        record = try_parse_line(with_field(6, "0M4I"))
        err = XAM2.cigar(record)
        @test err isa XAMError
        @test err.kind == XAMErrors.invalid_cigar
        @test err.wrapped_error == CIGARError(2, CIGARStrings.Errors.ZeroLength)
        @test mapping_span(record) == err
        @test mapping_blocks(record) == err
        @test occursin("in CIGAR", sprint(showerror, err))
        @test sprint(show, MIME"text/plain"(), record) isa String

        # So is the aux data
        record = try_parse_line(sam_line(VALID_FIELDS) * "\tbad")
        @test record isa SAMRecord
        @test aux_data(record) == XAMError(XAMErrors.invalid_aux_data)
        @test !XAM2.is_well_formed(aux_data_unchecked(record))
        @test occursin("Auxiliary: < malformed >", sprint(show, MIME"text/plain"(), record))
    end

    @testset "Invalid records" begin
        for (field, value, err) in [
                (1, "", XAMErrors.invalid_read_name),
                (1, "a b", XAMErrors.invalid_read_name),
                (1, "a@b", XAMErrors.invalid_read_name),
                (1, "!"^255, XAMErrors.invalid_read_name),
                (2, "x", XAMErrors.invalid_flags),
                (2, "-1", XAMErrors.invalid_flags),
                (2, "65536", XAMErrors.invalid_flags),
                (2, "99999999999999999999999", XAMErrors.invalid_flags),
                (3, "", XAMErrors.invalid_reference_name),
                (3, "*chr", XAMErrors.invalid_reference_name),
                (3, "=", XAMErrors.invalid_reference_name),
                (3, "a,b", XAMErrors.invalid_reference_name),
                (4, "-1", XAMErrors.invalid_position),
                (4, "1x", XAMErrors.invalid_position),
                (4, "2147483648", XAMErrors.invalid_position),
                (5, "256", XAMErrors.invalid_mapq),
                (5, "", XAMErrors.invalid_mapq),
                (6, "4", XAMErrors.invalid_cigar),
                (6, "M", XAMErrors.invalid_cigar),
                (6, "4Q", XAMErrors.invalid_cigar),
                (7, "==", XAMErrors.invalid_next_reference_name),
                (7, "*x", XAMErrors.invalid_next_reference_name),
                (8, "x", XAMErrors.invalid_next_position),
                (8, "2147483648", XAMErrors.invalid_next_position),
                (9, "--5", XAMErrors.invalid_template_length),
                (9, "2147483648", XAMErrors.invalid_template_length),
                (9, "-2147483649", XAMErrors.invalid_template_length),
                (9, "-", XAMErrors.invalid_template_length),
                (10, "AC1T", XAMErrors.invalid_sequence),
                (10, "", XAMErrors.invalid_sequence),
                (11, "II I", XAMErrors.invalid_quality),
                (11, "", XAMErrors.invalid_quality),
                # Structural checks of SEQ and QUAL
                (11, "III", XAMErrors.qual_seq_mismatch),
                (10, "*", XAMErrors.qual_no_seq),
            ]
            @test error_kind(try_parse_line(with_field(field, value))) == err
        end

        line = sam_line(VALID_FIELDS)
        @test error_kind(try_parse_line("")) == XAMErrors.missing_fields
        @test error_kind(try_parse_line(sam_line(VALID_FIELDS[1:10]))) == XAMErrors.missing_fields
        @test error_kind(try_parse_line(sam_line(VALID_FIELDS[1:3]))) == XAMErrors.missing_fields
        @test error_kind(try_parse_line(line * "\t")) == XAMErrors.invalid_aux_data
        @test error_kind(try_parse_line(line * "\tXA:Z:é")) == XAMErrors.invalid_aux_data
        @test error_kind(try_parse_line(line * "\r")) == XAMErrors.invalid_quality
        @test error_kind(try_parse_line(line * "\n")) == XAMErrors.invalid_quality
        @test error_kind(try_parse_line(replace(line, "\tACGT" => "\tAC\nGT"))) == XAMErrors.invalid_sequence

        # Columns: The invalid byte for syntax errors, the start of the field for
        # out-of-range values, and the end of the line if it ends early
        column(line) = try_parse_line(line).column
        # Column of the first byte of field `i`
        field_start(i) = ncodeunits(sam_line(VALID_FIELDS[1:(i - 1)])) + 2
        @test field_start(4) == 12
        @test column(with_field(4, "12x4")) == field_start(4) + 2
        @test column(with_field(4, "2147483648")) == field_start(4)
        @test column(with_field(1, "!"^255)) == 1
        @test column(with_field(11, "III")) == field_start(11)
        @test column(sam_line(VALID_FIELDS[1:3])) == ncodeunits(sam_line(VALID_FIELDS[1:3])) + 1
        @test try_parse_line(with_field(4, "-1")) == XAMError(XAMErrors.invalid_position; column = 12)
        @test sprint(showerror, try_parse_line(with_field(4, "-1"))) ==
            "XAMError: invalid position (POS) at column 12"

        # A failed parse resets the record, even if it was previously valid
        record = try_parse_line(line)
        @test error_kind(XAM2.try_parse!(record, ImmutableMemoryView(with_field(5, "256")))) == XAMErrors.invalid_mapq
        @test is_reset(record)
        @test_throws XAMError XAM2.parse!(record, ImmutableMemoryView(with_field(5, "256")))
    end

    @testset "Reusing records" begin
        # Regression: The cached CIGAR was not cleared when parsing a new record
        record = try_parse_line(sam_line(VALID_FIELDS))
        @test String(MemoryView(XAM2.cigar(record))) == "4M"
        XAM2.parse!(record, ImmutableMemoryView(with_field(6, "2M2S")))
        @test String(MemoryView(XAM2.cigar(record))) == "2M2S"

        # Records parsed without aux data after one with aux data
        XAM2.parse!(record, ImmutableMemoryView(sam_line(VALID_FIELDS) * "\tNM:i:1"))
        XAM2.parse!(record, ImmutableMemoryView(sam_line(VALID_FIELDS)))
        @test isempty(aux_data(record))
        @test quality_string(record) == "IIII"

        cp = copy(record)
        @test cp.data !== record.data
        @test XAM2.cigar(cp) == XAM2.cigar(record)
        @test query_name(cp) == "r1"
    end

    @testset "All errors can be shown" begin
        for kind in instances(XAMErrorKind)
            @test sprint(showerror, XAMError(kind)) isa String
        end
    end
end

summarize(r::SAMRecord) = (
    query_name(r),
    Integer(flags(r)),
    reference_name(r),
    left_position(r),
    next_reference_name(r),
    sequence(r) === nothing ? nothing : string(sequence(r)),
)

# A reader that reads `data`, and whose buffer holds at most `max_size` bytes
mutable struct LimitedBufReader <: XAM2.BufferIO.AbstractBufReader
    const data::Vector{UInt8}
    const max_size::Int
    # The buffer is `data[start:stop]`
    start::Int
    stop::Int
end

LimitedBufReader(s::String, max_size::Int) = LimitedBufReader(collect(codeunits(s)), max_size, 1, 0)

XAM2.BufferIO.get_buffer(x::LimitedBufReader) = ImmutableMemoryView(x.data)[x.start:x.stop]
XAM2.BufferIO.consume(x::LimitedBufReader, n::Int) = (x.start += n; nothing)
Base.close(::LimitedBufReader) = nothing

function XAM2.BufferIO.fill_buffer(x::LimitedBufReader)::Union{Int, Nothing}
    x.stop == length(x.data) && return 0
    n_buffered = x.stop - x.start + 1
    n_buffered == x.max_size && return nothing
    n = min(x.max_size - n_buffered, length(x.data) - x.stop)
    x.stop += n
    return n
end

@testset "SAMReader" begin
    bam_records = read_bam(EXAMPLE_BAM)
    records = SAMReader(collect, open(EXAMPLE_SAM))
    @test records isa Vector{SAMRecord}
    @test map(summarize, records) == map(summarize, bam_records)
    @test map(mapping_span, records) == map(mapping_span, bam_records)
    for f in (next_position, mapping_quality, template_length)
        @test map(f, records) == map(f, bam_records)
    end
    @test map(quality_string, records) ==
        [q === nothing ? nothing : String(q .+ UInt8(33)) for q in map(quality_scores, bam_records)]

    SAMReader(open(EXAMPLE_SAM)) do reader
        header = get_header(reader)
        @test References(header) == references(BAMReader(open(EXAMPLE_BAM); n_workers = 1))
        @test only(header.comments) == "An example comment"
        @test MemoryView(header) == MemoryView(get_header(BAMReader(open(EXAMPLE_BAM); n_workers = 1)))
    end

    # Copying BAM records to SAM records reproduces the SAM lines
    sam_lines = filter(l -> !startswith(l, '@'), readlines(EXAMPLE_SAM))
    record = SAMRecord()
    @test [String(copy(copy!(record, r).data)) for r in bam_records] == sam_lines

    # CRLF line endings, and no trailing newline
    data = "@HD\tVN:1.6\r\n" * join(sam_lines, "\r\n")
    @test map(summarize, SAMReader(collect, IOBuffer(data))) == map(summarize, records)

    # Header only, and empty file
    @test isempty(SAMReader(collect, IOBuffer("@HD\tVN:1.6\n")))
    @test isempty(SAMReader(collect, IOBuffer("")))

    # Invalid header
    @test_throws XAMError SAMReader(IOBuffer("@HD\tVN:1.6\n@XX\n"))

    @testset "Errors are sticky, and the record stays valid" begin
        data = "@HD\tVN:1.6\n@CO\tx\n" *
            join([sam_line(VALID_FIELDS), with_field(4, "-1"), sam_line(VALID_FIELDS)], '\n')
        expected = XAMError(XAMErrors.invalid_position; record = 2, line = 4, column = 12)
        path = tempname()
        write(path, data)
        reader = SAMReader(open(path); copy = false)
        record = XAM2.try_read!(reader)
        @test query_name(record) == "r1"
        @test XAM2.try_read!(reader) == expected
        @test is_reset(record)
        @test XAM2.try_read!(reader) === expected
        @test sprint(showerror, expected) ==
            "XAMError: invalid position (POS) at record 2, line 4, column 12"
        @test_throws XAMError iterate(reader)

        # Closing sets a permanent error
        close(reader)
        @test error_kind(XAM2.try_read!(reader)) == XAMErrors.closed_stream
    end

    @testset "try_init" begin
        reader = XAM2.try_init(SAMReader, open(EXAMPLE_SAM))
        @test reader isa SAMReader
        @test length(collect(reader)) == 6
        close(reader)

        # On failure, the error is returned, and the input is closed
        io = IOBuffer("@HD\tVN:1.6\n@SQ\tSN:x\n" * sam_line(VALID_FIELDS))
        @test XAM2.try_init(SAMReader, io) == XAMError(XAMErrors.header_missing_ln; line = 2, column = 1)
        @test !isopen(io)
        io = IOBuffer("@HD\tVN:1.6\n@SQ\tSN:x\n")
        @test_throws XAMError(XAMErrors.header_missing_ln; line = 2, column = 1) SAMReader(io)
        @test !isopen(io)
    end

    @testset "Readers that can't buffer a whole line" begin
        header = "@HD\tVN:1.6\n@CO\t" * "x"^100 * "\n"
        short = sam_line(VALID_FIELDS)
        long = sam_line(["r2"; VALID_FIELDS[2:9]; "A"^100; "I"^100])
        data = header * join([short, long, short], '\n')
        expected = map(summarize, SAMReader(collect, IOBuffer(data)))
        @test length(expected) == 3

        # Lines longer than the buffer, including header lines, and CRLF line endings,
        # where the buffer may end between the CR and LF.
        for text in (data, replace(data, '\n' => "\r\n"))
            for max_size in (1, 2, 7, length(short), length(short) + 1, length(long), length(long) + 1)
                SAMReader(LimitedBufReader(text, max_size)) do reader
                    @test only(get_header(reader).comments) == "x"^100
                    @test map(summarize, reader) == expected
                end
            end
        end

        @testset "Maximum line length" begin
            get_line(s, max_length) = XAM2.get_next_line_view(XAM2.BufferIO.ExpandReader(LimitedBufReader(s, 2)), max_length)
            @test get_line("abcd\nef", 4) == ImmutableMemoryView("abcd\n")
            @test get_line("abcd\r\nef", 4) == ImmutableMemoryView("abcd\r\n")
            @test get_line("abcd", 4) == ImmutableMemoryView("abcd")
            @test get_line("", 4) === nothing
            # Without newlines, buffering stops when the line is too long
            @test get_line("a"^10_000, 4) == XAMErrors.line_too_long
            @test get_line("a"^10_000 * "\n", 100) == XAMErrors.line_too_long
        end

        @testset "Bgzipped SAM" begin
            mktempdir() do dir
                # Large enough that lines cross BGZF block boundaries
                path = joinpath(dir, "file.sam.gz")
                n_records = 10_000
                BGZFWriter(open(path, "w")) do writer
                    write(writer, header)
                    for _ in 1:n_records
                        write(writer, short, '\n')
                    end
                end
                records = SAMReader(collect, BGZFReader(open(path); n_workers = 1))
                @test length(records) == n_records
                @test all(r -> summarize(r) == expected[1], records)

                # Block boundaries inside a header line and a record line.
                # Each chunk is compressed to its own blocks.
                function bgzf_chunks(splits...)
                    bounds = [0; collect(splits); ncodeunits(data)]
                    chunks = [codeunits(data)[(bounds[i] + 1):bounds[i + 1]] for i in 1:(length(bounds) - 1)]
                    return [bgzf_compress(dir, collect(chunk)) for chunk in chunks]
                end
                record_start = ncodeunits(header)
                bytes = reduce(vcat, bgzf_chunks(5, record_start + 10))
                SAMReader(BGZFReader(IOBuffer(bytes); n_workers = 1)) do reader
                    @test only(get_header(reader).comments) == "x"^100
                    @test map(summarize, reader) == expected
                end

                # Malformed data in the block after the line's first block is thrown
                long_start = record_start + ncodeunits(short) + 1
                blocks = bgzf_chunks(long_start + 10)
                corrupt_crc!(blocks[2])
                reader = SAMReader(BGZFReader(IOBuffer(reduce(vcat, blocks)); n_workers = 1))
                @test query_name(XAM2.try_read!(reader)) == "r1"
                @test_throws BGZFError XAM2.try_read!(reader)
                close(reader)
            end
        end
    end
end
