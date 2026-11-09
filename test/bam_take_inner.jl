@testset "BAM take_inner" begin
    (header, records) = BAMReader(open(EXAMPLE_BAM); n_workers = 1) do reader
        (get_header(reader)::SAMHeader, collect(reader))
    end

    @testset "BAMWriter" begin
        # Writing to an IOBuffer, which would be lost by closing it
        buffer = IOBuffer()
        bgzf = BAMWriter(buffer, header; n_workers = 1) do writer
            foreach(r -> write(writer, r), records)
            bgzf = @inferred take_inner(writer)
            @test !isopen(writer)
            @test_throws IOError write(writer, first(records))
            @test_throws IOError flush(writer)
            @test_throws IOError take_inner(writer)
            bgzf
        end
        # Closing after `take_inner` does not close the BGZFWriter
        @test bgzf isa BGZFWriter{BufWriter{IOBuffer}}
        @test isopen(bgzf)
        # Taking the BufWriter from the BGZFWriter writes the EOF block,
        # so the output is a complete BAM file
        flush(take_inner(bgzf))
        @test BAMReader(collect, IOBuffer(take!(buffer)); n_workers = 1) == records

        # From a BGZFWriter, the same BGZFWriter is returned
        bgzf = BGZFWriter(IOBuffer(); n_workers = 1)
        @test take_inner(BAMWriter(bgzf, header)) === bgzf
        close(bgzf)
    end

    # The records read after seeking a new reader to `offset`
    records_from(offset) = BAMReader(open(EXAMPLE_BAM); n_workers = 1) do reader
        collect(seek(reader, offset))
    end

    @testset "BAMReader" begin
        io = open(EXAMPLE_BAM)
        bgzf = BAMReader(io; n_workers = 2, copy = false) do reader
            record = XAM2.try_read!(reader)
            offset = virtual_position(reader)
            bgzf = @inferred take_inner(reader)
            @test !isopen(reader)
            @test XAM2.try_read!(reader) == XAMError(XAMErrors.closed_stream)
            @test_throws IOError take_inner(reader)
            @test_throws IOError seek_first_record(reader)
            @test_throws IOError virtual_position(reader)
            # The record is still valid
            @test record == first(records)
            # Positioned after the first record
            @test virtual_position(bgzf) == offset
            @test records_from(offset) == records[2:end]
            bgzf
        end
        # Closing after `take_inner` does not close the BGZFReader
        @test bgzf isa BGZFReader{BufReader{IOStream}}
        @test isopen(bgzf)
        @test isopen(io)
        close(bgzf)
        @test !isopen(io)

        # From a BGZFReader, the same BGZFReader is returned, positioned at the
        # first record, and EOF after all records
        bgzf = BGZFReader(open(EXAMPLE_BAM); n_workers = 1)
        reader = BAMReader(bgzf)
        offset = virtual_position(reader)
        @test take_inner(reader) === bgzf
        @test virtual_position(bgzf) == offset
        @test records_from(offset) == records
        close(bgzf)
        reader = BAMReader(open(EXAMPLE_BAM); n_workers = 1)
        collect(reader)
        bgzf = take_inner(reader)
        @test eof(bgzf)
        close(bgzf)
    end
end
