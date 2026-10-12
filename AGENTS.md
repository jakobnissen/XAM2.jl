This repo, XAM2.jl is an experimental rewrite of the package XAM.jl, a reader/writer of bioinformatics SAM and BAM files.
XAM2's rewrite is motivated by several factors which inform its design space:
* The packages BufferIO, MemoryViews and BGZFLib, which are useful for XAM, and higher performance than previous alternatives
* Poor adherence to the SAM spec in the existing XAM package
* Subpar performance, where a XAM library ought to have excellent performance.

Currently, the code is in a work-in-progress state, and is being dogfooded by an application called Overburden.jl.
Rapid and breaking changes are expected and need not be documented.
Only easy and superficial tests should be added; it is expected the code will be refactored anyway.

## Design criteria
* Use the type system to cover spec edge cases, e.g. return `nothing` for missing fields.
* Distinguish between in-data edge cases (like missing data) versus IO or filesystem errors. The former should be represented
  in the file system, the latter as exceptions, as they can rarely be sensibly handled.
* Performance is very important for this codebase.

## Coding guidelines
After making changes, format with `runic -i .`. If executable is not found or errors, alert the user, and do not attempt to install Runic.
For all user-facing operations, there must be a version that does not throw on malformed data. Bad user-provided inputs may cause exceptions.
Records are parsed semi-lazily: Simple structural invariants are checked on record construction, which allows `@inbounds` to be used safely. However, variable-length data is lazily checked.
