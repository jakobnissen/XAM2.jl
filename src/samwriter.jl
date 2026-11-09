"""
    SAMWriter(io::Union{IO, AbstractBufWriter}, header::SAMHeader)::SAMWriter

Construct a `SAMWriter`, which writes `SAMRecord`s writing to `io`.
If `io` is an `IO`, it is wrapped in a `BufferIO.BufWriter`.
The text of `header` is written on instantiation, followed by a newline if
it does not end with one.

The writer takes ownership of `io`, so `io` is not to be mutated
while the writer is used.
Closing the resulting writer will close `io`.

See [`XAM2.try_write`](@ref).

# Examples
```julia
BAMReader(open("in.bam")) do reader
    header = get_header(reader)
    header isa XAMError && throw(header)
    SAMWriter(open("out.sam", "w"), header) do writer
        for record in reader
            is_mapped(record) && write(writer, record)
        end
    end
end
```

# Extended help
The reference names of records are not checked against
the `@SQ` lines of `header`.

"""
mutable struct SAMWriter{T <: AbstractBufWriter}
    const writer::T
    # Scratch record used to convert `BAMRecord`s to SAM
    const record::SAMRecord
    # Set by `close` and `take_inner`
    closed::Bool
end

function SAMWriter(io::IO, header::SAMHeader)
    return SAMWriter(BufWriter(io), header)
end

function SAMWriter(writer::T, header::SAMHeader) where {T <: AbstractBufWriter}
    text = MemoryView(header)
    write(writer, text)
    if !isempty(text) && @inbounds(text[end]) != UInt8('\n')
        write(writer, UInt8('\n'))
    end
    return SAMWriter{T}(writer, SAMRecord(), false)
end

function SAMWriter(f::Function, args::Vararg{Any, N}; kwargs...) where {N}
    writer = SAMWriter(args...; kwargs...)
    return try
        f(writer)
    finally
        close(writer)
    end
end

function Base.show(io::IO, writer::SAMWriter)
    print(io, '<')
    summary(io, writer)
    return print(io, '>')
end

# Closing twice is a no-op, so that closing after `take_inner` does not close the
# taken writer, e.g. in the `do`-block form of `SAMWriter`.
function Base.close(writer::SAMWriter)
    writer.closed && return nothing
    writer.closed = true
    return close(writer.writer)
end

Base.flush(writer::SAMWriter) = flush(writer.writer)
Base.isopen(writer::SAMWriter) = !writer.closed && isopen(writer.writer)

"""
    BufferIO.take_inner(writer::SAMWriter{T})::T

Close `writer` without closing or flushing its underlying writer, and return it.
If `writer` was constructed from an `IO`, this is the `BufferIO.BufWriter` wrapping it.
Since data is written directly to the underlying writer, no written data is lost.

Calling this on a closed `writer` throws a `BufferIO.IOError` of kind `ClosedIO`.
After this, writing to `writer` throws the same error, and `close(writer)` does nothing.

# Examples
```julia
# Write SAM to stdout without closing it
writer = SAMWriter(stdout, header)
write(writer, record)
flush(BufferIO.take_inner(writer))
```
"""
function BufferIO.take_inner(writer::SAMWriter{T})::T where {T}
    writer.closed && throw(IOError(IOErrorKinds.ClosedIO))
    writer.closed = true
    return writer.writer
end

"""
    XAM2.try_write(writer::SAMWriter, record::Union{SAMRecord, BAMRecord})::Union{Int, XAMError}

Write `record` as a line to `writer`, and return the number of bytes written.
For `SAMRecord`s, its `data` field is written without validation.
A `BAMRecord` is first converted to SAM. On conversion error,
write nothing and return the `XAMError`,

See also: [`SAMWriter`](@ref), [`XAM2.try_copy!`](@ref)
"""
function try_write(writer::SAMWriter, record::SAMRecord)::Union{Int, XAMError}
    return write_line(writer, record)
end

function try_write(writer::SAMWriter, record::BAMRecord)::Union{Int, XAMError}
    sam = try_copy!(writer.record, record)
    sam isa XAMError && return sam
    return write_line(writer, sam)
end

function write_line(writer::SAMWriter, record::SAMRecord)::Int
    writer.closed && throw(IOError(IOErrorKinds.ClosedIO))
    io = writer.writer
    result = write(io, record.data)
    return result + write(io, UInt8('\n'))
end

function Base.write(writer::SAMWriter, record::Union{SAMRecord, BAMRecord})
    result = try_write(writer, record)
    result isa XAMError && throw(result)
    return result
end
