# Trim smoke test: A program exercising XAM2's public API, compiled with
# `juliac --trim=safe` by build.sh. It reads a BAM file twice, writes it with
# the references in reverse order (exercising reference ID translation), and
# checks the result. It then reads the same records from a SAM file.
# Exits with 0 on success.

using XAM2:
    XAM2,
    BAMReader,
    BAMWriter,
    SAMReader,
    References,
    SAMHeader,
    XAMErrorKind,
    references,
    query_name,
    reference_name,
    next_reference_name,
    mapping_span,
    try_write

function summarize(reader::Union{BAMReader, SAMReader})
    n_records = 0
    names = String[]
    for record in reader
        n_records += 1
        name = query_name(record)
        rname = reference_name(record)
        rnext = next_reference_name(record)
        push!(names, name === nothing ? "*" : String(name))
        push!(names, rname === nothing ? "*" : String(rname))
        push!(names, rnext === nothing ? "*" : String(rnext))
        span = mapping_span(record)
        push!(names, span isa UnitRange{Int} ? string(span) : "*")
    end
    return (n_records, names)
end

function (@main)(args::Vector{String})::Cint
    (in_path, out_path, sam_path) = (args[1], args[2], args[3])

    reader = BAMReader(open(in_path); n_workers = 1)
    (n1, names1) = summarize(reader)
    XAM2.seek_first_record(reader)
    (n2, names2) = summarize(reader)

    # Write all records, with the references in reverse order
    refs = references(reader)
    reversed = References(
        String[String(ref.name) for ref in Iterators.reverse(refs)],
        Int[ref.length for ref in Iterators.reverse(refs)],
    )
    header = XAM2.try_parse_header(XAM2.ImmutableMemoryView("@HD\tVN:1.6\n"))
    header isa SAMHeader || return 1
    XAM2.seek_first_record(reader)
    writer = BAMWriter(open(out_path, "w"), header, reversed; n_workers = 1)
    n_errors = 0
    for record in reader
        n_errors += try_write(writer, record) isa XAMErrorKind
    end
    close(writer)
    close(reader)

    reader = BAMReader(open(out_path); n_workers = 1)
    (n3, names3) = summarize(reader)
    close(reader)

    reader = SAMReader(open(sam_path))
    (n4, names4) = summarize(reader)
    close(reader)

    println(Core.stdout, "records: ", n1, " ", n2, " ", n3, " ", n4)
    println(Core.stdout, "write errors: ", n_errors)
    ok = n1 > 0 && n1 == n2 == n3 == n4 && iszero(n_errors) && names1 == names2 == names3 == names4
    println(Core.stdout, ok ? "OK" : "FAILED")
    return ok ? 0 : 1
end
