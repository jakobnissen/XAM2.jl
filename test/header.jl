# Tests of parsing and validating SAM headers


# Return (line, column, kind) of the error from parsing `s`
function header_error(s::String)
    result = try_parse_header(ImmutableMemoryView(s))
    result isa XAMError || return result
    return (result.line, result.column, result.kind)
end

@testset "SAM header" begin
    @testset "Valid header" begin
        header = parse_header(
            "@HD\tVN:1.6\tSO:coordinate\tGO:none\n" *
                "@SQ\tSN:chr1\tLN:1000\tTP:circular\tAS:hg38\n" *
                "@SQ\tSN:HLA-A*01:01\tLN:2147483647\n" *
                "@RG\tID:rg1\tSM:sample 1\n" *
                "@PG\tID:prog\tPN:prog\tCL:prog --flag\n" *
                "@CO\tA comment\twith tabs: Æxámplè\n" *
                "@CO\n" *
                "@RG\tID:rg1\n" # duplicate IDs are not validated
        )
        hd = header.header
        @test (hd.major, hd.minor) == (1, 6)
        @test hd.sort_order == XAM2.SortOrders.coordinate
        @test hd["GO"] == "none"
        @test length(hd) == 3

        (sq1, sq2) = header.seqs
        @test sq1.sn == "chr1"
        @test sq1.sn isa AbstractString
        @test sq1.ln === 1000
        @test sq1.topology == XAM2.Topologies.circular
        @test sq1["AS"] == "hg38"
        @test collect(sq1) == [
            AuxTag("SN") => "chr1", AuxTag("LN") => "1000", AuxTag("TP") => "circular", AuxTag("AS") => "hg38",
        ]
        @test get(sq1, "XX", nothing) === nothing
        @test String(MemoryView(sq1)) == "SN:chr1\tLN:1000\tTP:circular\tAS:hg38"
        @test sq2.sn == "HLA-A*01:01"
        @test sq2.ln == typemax(Int32)
        @test sq2.topology == XAM2.Topologies.linear

        @test [rg.id for rg in header.read_groups] == ["rg1", "rg1"]
        @test header.read_groups[1]["SM"] == "sample 1"
        @test only(header.programs).id == "prog"
        @test header.programs[1]["CL"] == "prog --flag"
        @test header.comments == ["A comment\twith tabs: Æxámplè", ""]
        @test all(c -> c isa AbstractString, header.comments)

        # Defaults
        header = parse_header("@SQ\tSN:x\tLN:1")
        @test header.header === nothing
        @test only(header.seqs).sn == "x"
        @test parse_header("@HD\tVN:1.6\n").header.sort_order == XAM2.SortOrders.unknown

        # CRLF line endings, and no trailing newline
        header = parse_header("@HD\tVN:1.6\r\n@SQ\tSN:x\tLN:10\r\n@CO\tc\r\n@RG\tID:a\r")
        @test only(header.seqs).ln == 10
        @test only(header.comments) == "c"
        @test only(header.read_groups).id == "a"

        # Empty header
        header = parse_header("")
        @test header.header === nothing
        @test isempty(header.seqs)
        @test isempty(MemoryView(header))
    end

    @testset "Invalid header" begin
        E = XAMErrors
        for (s, err) in [
                # Line structure
                "SQ\tSN:x\tLN:1" => (1, 1, E.header_invalid_line),
                "@SQ SN:x\tLN:1" => (1, 1, E.header_invalid_line),
                "@S" => (1, 1, E.header_invalid_line),
                "@HD\tVN:1.6\n\n@SQ\tSN:x\tLN:1" => (2, 1, E.header_invalid_line),
                "\n" => (1, 1, E.header_invalid_line),
                "\r\n" => (1, 1, E.header_invalid_line),
                "@HD\tVN:1.6\n\n" => (2, 1, E.header_invalid_line),
                "@XY\tID:1" => (1, 2, E.header_unknown_record_type),
                "@hd\tVN:1.6" => (1, 2, E.header_unknown_record_type),
                "@SQ\tSN:x\t\tLN:1" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\tLN:1\t" => (1, 15, E.header_invalid_field),
                "@SQ\t" => (1, 5, E.header_invalid_field),
                "@SQ\tSN:x\tLN" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\tLN:" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\t1N:1" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\tL_:1" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\tLN=1" => (1, 10, E.header_invalid_field),
                "@SQ\tSN:x\tLN:1\tSN:y" => (1, 15, E.header_duplicate_tag),
                "@RG\tID:a\tSM:b\tSM:b" => (1, 15, E.header_duplicate_tag),
                # HD
                "@SQ\tSN:x\tLN:1\n@HD\tVN:1.6" => (2, 1, E.header_hd_not_first),
                "@HD\tVN:1.6\n@HD\tVN:1.6" => (2, 1, E.header_hd_not_first),
                "@HD\tSO:coordinate" => (1, 1, E.header_missing_vn),
                "@HD" => (1, 1, E.header_missing_vn),
                "@HD\tVN:1" => (1, 5, E.header_invalid_vn),
                "@HD\tVN:1.x" => (1, 5, E.header_invalid_vn),
                "@HD\tVN:.6" => (1, 5, E.header_invalid_vn),
                "@HD\tVN:1.6\tSO:sorted" => (1, 12, E.header_invalid_sort_order),
                # SQ
                "@SQ\tLN:1" => (1, 1, E.header_missing_sn),
                "@SQ\tSN:x" => (1, 1, E.header_missing_ln),
                "@SQ\tSN:*x\tLN:1" => (1, 5, E.header_invalid_sn),
                "@SQ\tSN:=x\tLN:1" => (1, 5, E.header_invalid_sn),
                "@SQ\tSN:a b\tLN:1" => (1, 5, E.header_invalid_sn),
                "@SQ\tSN:a,b\tLN:1" => (1, 5, E.header_invalid_sn),
                "@SQ\tSN:x\tLN:1\n@SQ\tSN:y\tLN:1\n@SQ\tLN:2\tSN:x" => (3, 10, E.header_duplicate_sn),
                "@SQ\tSN:x\tLN:0" => (1, 10, E.header_invalid_ln),
                "@SQ\tSN:x\tLN:-1" => (1, 10, E.header_invalid_ln),
                "@SQ\tSN:x\tLN:+1" => (1, 10, E.header_invalid_ln),
                "@SQ\tSN:x\tLN:2147483648" => (1, 10, E.header_invalid_ln),
                "@SQ\tSN:x\tLN:10x" => (1, 10, E.header_invalid_ln),
                "@SQ\tSN:x\tLN:1\tTP:round" => (1, 15, E.header_invalid_topology),
                # RG and PG
                "@RG\tSM:x" => (1, 1, E.header_missing_id),
                "@PG\tPN:x" => (1, 1, E.header_missing_id),
            ]
            @test header_error(s) == err
        end

        @test sprint(showerror, try_parse_header(ImmutableMemoryView("@HD\tVN:1.6\n@PG\tPN:x"))) ==
            "XAMError: @RG or @PG line has no ID field at line 2, column 1"
    end
end

header_text(header::SAMHeader) = String(MemoryView(header))

@testset "Editing headers" begin
    SO = XAM2.SortOrders
    @testset "add_program" begin
        add(s, id; kw...) = header_text(XAM2.add_program(parse_header(s), id; kw...))

        # Placed after the last @PG line, or where @PG lines conventionally go
        @test add("", "a") == "@PG\tID:a\n"
        @test add("@HD\tVN:1.6\n@SQ\tSN:x\tLN:1\n@CO\tc", "a"; pn = "A", cl = "a --x", ds = "d", vn = "1") ==
            "@HD\tVN:1.6\n@SQ\tSN:x\tLN:1\n@PG\tID:a\tPN:A\tCL:a --x\tDS:d\tVN:1\n@CO\tc\n"
        @test add("@PG\tID:x\n@CO\tc\n@PG\tID:y\tPP:x\n@CO\td", "a") ==
            "@PG\tID:x\n@CO\tc\n@PG\tID:y\tPP:x\n@PG\tID:a\tPP:y\n@CO\td\n"
        @test add("@CO\tc\n@RG\tID:r", "a") == "@CO\tc\n@RG\tID:r\n@PG\tID:a\n"

        # CRLF line endings are normalized, and the original header is unchanged
        header = parse_header("@HD\tVN:1.6\r\n@PG\tID:x")
        @test header_text(XAM2.add_program(header, "a")) == "@HD\tVN:1.6\n@PG\tID:x\n@PG\tID:a\tPP:x\n"
        @test header_text(header) == "@HD\tVN:1.6\r\n@PG\tID:x"

        # Unique IDs
        @test add("@PG\tID:a", "a") == "@PG\tID:a\n@PG\tID:a.1\tPP:a\n"
        @test add("@PG\tID:a\n@PG\tID:a.1\tPP:a", "a"; pp = nothing) ==
            "@PG\tID:a\n@PG\tID:a.1\tPP:a\n@PG\tID:a.2\n"
        @test last(XAM2.add_program(parse_header("@PG\tID:a"), "a").programs).id == "a.1"

        # Default PP is the end of the chain, if there is exactly one chain
        @test add("@PG\tID:x\tPP:y\n@PG\tID:y", "a") == "@PG\tID:x\tPP:y\n@PG\tID:y\n@PG\tID:a\tPP:x\n"
        @test add("@PG\tID:x\n@PG\tID:y", "a") == "@PG\tID:x\n@PG\tID:y\n@PG\tID:a\n"
        @test add("@PG\tID:x\tPP:y\n@PG\tID:y\tPP:x", "a") == "@PG\tID:x\tPP:y\n@PG\tID:y\tPP:x\n@PG\tID:a\n"
        @test add("@PG\tID:x\tPP:nonexistent", "a") == "@PG\tID:x\tPP:nonexistent\n@PG\tID:a\tPP:x\n"
        # Duplicate IDs in the input count as one chain end
        @test add("@PG\tID:x\n@PG\tID:x", "a") == "@PG\tID:x\n@PG\tID:x\n@PG\tID:a\tPP:x\n"

        # Explicit PP
        @test add("@PG\tID:x\n@PG\tID:y", "a"; pp = "x") == "@PG\tID:x\n@PG\tID:y\n@PG\tID:a\tPP:x\n"
        @test add("@PG\tID:x", "a"; pp = nothing) == "@PG\tID:x\n@PG\tID:a\n"
        @test_throws ArgumentError add("@PG\tID:x", "a"; pp = "y")
        @test_throws ArgumentError add("", "a"; pp = "a")

        # Invalid values
        for bad in ["", "a\tb", "a\nb", "a\rb", "a\0b"]
            @test_throws ArgumentError add("", bad)
            @test_throws ArgumentError add("", "a"; cl = bad)
        end
    end

    @testset "add_read_group" begin
        add(s, id; kw...) = header_text(XAM2.add_read_group(parse_header(s), id; kw...))
        @test add("", "r") == "@RG\tID:r\n"
        @test add("@HD\tVN:1.6\n@SQ\tSN:x\tLN:1\n@PG\tID:p\n@CO\tc", "r"; sm = "s", lb = "l", pi = "300") ==
            "@HD\tVN:1.6\n@SQ\tSN:x\tLN:1\n@RG\tID:r\tLB:l\tPI:300\tSM:s\n@PG\tID:p\n@CO\tc\n"
        @test add("@RG\tID:a\n@PG\tID:p\n@RG\tID:b", "r") == "@RG\tID:a\n@PG\tID:p\n@RG\tID:b\n@RG\tID:r\n"
        header = XAM2.add_read_group(parse_header(""), "r"; sm = "sample 1", pl = "ILLUMINA")
        @test only(header.read_groups)["SM"] == "sample 1"

        @test_throws ArgumentError add("@RG\tID:r", "r")
        for bad in ["", "a\tb", "a\nb", "a\rb", "a\0b"]
            @test_throws ArgumentError add("", bad)
            @test_throws ArgumentError add("", "r"; sm = bad)
        end
    end

    @testset "set_sort_order" begin
        set(s, so) = header_text(XAM2.set_sort_order(parse_header(s), so))
        @test set("", SO.coordinate) == "@HD\tVN:1.6\tSO:coordinate\n"
        @test set("@SQ\tSN:x\tLN:1", SO.queryname) == "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:x\tLN:1\n"
        @test set("@HD\tVN:1.5\n@CO\tc", SO.unsorted) == "@HD\tVN:1.5\tSO:unsorted\n@CO\tc\n"
        # SO is replaced in place, and SS and GO are removed
        @test set("@HD\tSO:queryname\tVN:1.6\tSS:queryname:natural\tGO:query\tXY:z", SO.unknown) ==
            "@HD\tSO:unknown\tVN:1.6\tXY:z\n"
        for so in instances(SO.SortOrder)
            @test XAM2.set_sort_order(parse_header(""), so).header.sort_order == so
        end
    end

    @testset "with_references" begin
        with(s, names, lengths) = header_text(XAM2.with_references(parse_header(s), References(names, lengths)))
        # Lines with same name and length are kept, others replaced
        @test with("@HD\tVN:1.6\n@SQ\tSN:a\tLN:10\tM5:x\n@SQ\tLN:20\tSN:b\tAS:y\n@SQ\tSN:c\tLN:1\n@RG\tID:r", ["c", "b", "a", "d"], [1, 20, 11, 5]) ==
            "@HD\tVN:1.6\n@SQ\tSN:c\tLN:1\n@SQ\tLN:20\tSN:b\tAS:y\n@SQ\tSN:a\tLN:11\n@SQ\tSN:d\tLN:5\n@RG\tID:r\n"
        # New lines are placed at the first @SQ line, or after the @HD line
        @test with("@RG\tID:r\n@SQ\tSN:a\tLN:1\n@CO\tc\n@SQ\tSN:b\tLN:1", ["b"], [1]) ==
            "@RG\tID:r\n@SQ\tSN:b\tLN:1\n@CO\tc\n"
        @test with("@HD\tVN:1.6\n@RG\tID:r", ["a"], [1]) == "@HD\tVN:1.6\n@SQ\tSN:a\tLN:1\n@RG\tID:r\n"
        @test with("@RG\tID:r", ["a"], [1]) == "@SQ\tSN:a\tLN:1\n@RG\tID:r\n"
        @test with("@HD\tVN:1.6\n@SQ\tSN:a\tLN:1", String[], Int[]) == "@HD\tVN:1.6\n"

        refs = References(["x", "y"], [5, 6])
        @test References(XAM2.with_references(parse_header("@SQ\tSN:a\tLN:1"), refs)) == refs

        @test_throws ArgumentError with("", ["a b"], [1])
        @test_throws ArgumentError with("", ["*a"], [1])
    end

    @testset "Writing edited header" begin
        reader = BAMReader(open(EXAMPLE_BAM); n_workers = 1)
        header = get_header(reader)::SAMHeader
        close(reader)
        refs = References(["chr2", "chrX"], [500, 100])
        header = XAM2.with_references(header, refs)
        header = XAM2.set_sort_order(header, SO.unsorted)
        header = XAM2.add_read_group(header, "rg2"; sm = "sample2")
        header = XAM2.add_program(header, "test"; pn = "test", cl = "test --flag")
        mktempdir() do dir
            path = write_bam(dir, BAMRecord[], header, refs)
            reader = BAMReader(open(path); n_workers = 1)
            @test MemoryView(get_header(reader)) == MemoryView(header)
            @test references(reader) == refs
            close(reader)
        end
        @test header_text(header) ==
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr2\tLN:500\n@SQ\tSN:chrX\tLN:100\n" *
            "@RG\tID:rg1\tSM:sample1\n@RG\tID:rg2\tSM:sample2\n" *
            "@PG\tID:prog\tPN:prog\n@PG\tID:test\tPN:test\tCL:test --flag\tPP:prog\n" *
            "@CO\tAn example comment\n"
    end
end
