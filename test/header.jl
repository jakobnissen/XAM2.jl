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
