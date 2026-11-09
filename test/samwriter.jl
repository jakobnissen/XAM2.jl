# Write `records` with a `SAMWriter`, and return the written text
function write_sam(header::SAMHeader, records)
    io = IOBuffer()
    writer = SAMWriter(io, header)
    for record in records
        write(writer, record)
    end
    flush(take_inner(writer))
    return String(take!(io))
end

@testset "SAMWriter" begin
    sam_text = read(EXAMPLE_SAM, String)

    @testset "SAM to SAM is verbatim" begin
        (header, records) = SAMReader(open(EXAMPLE_SAM)) do reader
            (get_header(reader), collect(reader))
        end
        @test write_sam(header, records) == sam_text
    end

    @testset "BAM to SAM" begin
        (header, records) = BAMReader(open(EXAMPLE_BAM); n_workers = 1) do reader
            (get_header(reader)::SAMHeader, collect(reader))
        end
        @test write_sam(header, records) == sam_text
    end

    @testset "Header" begin
        # A newline is added to a header without one, and an empty header writes nothing
        @test write_sam(parse_header("@HD\tVN:1.6"), []) == "@HD\tVN:1.6\n"
        @test write_sam(parse_header("@HD\tVN:1.6\n"), []) == "@HD\tVN:1.6\n"
        @test write_sam(parse_header(""), [try_parse_line(sam_line(VALID_FIELDS))]) ==
            sam_line(VALID_FIELDS) * '\n'
    end

    @testset "Return values and errors" begin
        io = IOBuffer()
        writer = SAMWriter(io, parse_header(""))
        line = sam_line(VALID_FIELDS)
        @test try_write(writer, try_parse_line(line)) == ncodeunits(line) + 1

        # A BAM record that can't be converted writes nothing
        record = copy!(BAMRecord(References(["chr1", "chr2"], [1000, 500])), try_parse_line(line))
        push!(record.data, UInt8('X'))
        @test try_write(writer, record) == XAMError(XAMErrors.invalid_aux_data)
        @test_throws XAMError write(writer, record)
        flush(writer)
        @test String(take!(io)) == line * '\n'

        @test isopen(writer)
        close(writer)
        @test !isopen(writer)
        @test !isopen(io)
    end

    @testset "take_inner" begin
        io = IOBuffer()
        line = sam_line(VALID_FIELDS)
        inner = SAMWriter(io, parse_header("")) do writer
            write(writer, try_parse_line(line))
            inner = take_inner(writer)
            @test inner isa BufWriter
            @test !isopen(writer)
            @test_throws IOError write(writer, try_parse_line(line))
            @test_throws IOError take_inner(writer)
            inner
        end
        # Closing the writer after `take_inner` does not close the inner writer
        @test isopen(inner)
        write(inner, "more")
        flush(inner)
        @test String(take!(io)) == line * "\nmore"
    end
end

@testset "SAMReader take_inner" begin
    header = "@HD\tVN:1.6\n"
    line = sam_line(VALID_FIELDS)
    data = header * line * '\n' * line * "\nrest"
    bufreader = BufReader(IOBuffer(data))
    expander = ExpandReader(BufReader(IOBuffer(data)))
    # (input, expected type of the taken reader, whether it is the input)
    cases = [
        (IOBuffer(data), BufReader{IOBuffer}, false),
        (bufreader, BufReader{IOBuffer}, true),
        (expander, typeof(expander), true),
        # Not expandable, so it is returned wrapped in an `ExpandReader`
        (LimitedBufReader(data, 7), ExpandReader{LimitedBufReader}, false),
    ]
    for (input, T, is_input) in cases
        (record, inner) = SAMReader(input; copy = false) do reader
            record = XAM2.try_read!(reader)
            inner = @inferred take_inner(reader)
            @test !isopen(reader)
            @test XAM2.try_read!(reader) == XAMError(XAMErrors.closed_stream)
            @test_throws IOError take_inner(reader)
            (record, inner)
        end
        @test inner isa T
        is_input && @test inner === input
        # The inner reader is not closed, keeps its buffered data,
        # and is positioned after the last record read
        @test isopen(inner)
        @test read(inner, String) == line * "\nrest"
        # The reader's record is still valid
        @test query_name(record) == "r1"
    end

    # After the header, if no records were read
    reader = SAMReader(IOBuffer(data))
    @test read(take_inner(reader), String) == line * '\n' * line * "\nrest"
end
