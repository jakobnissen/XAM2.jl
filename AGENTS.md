After making changes, format with `runic -i .`. If executable is not found or errors, alert the user, and do not attempt to install Runic.
For all user-facing operations, there must be a version that does not throw on malformed data. Bad user-provided inputs may cause exceptions.
Records are parsed semi-lazily: Simple structural invariants are checked on record construction, which allows `@inbounds` to be used safely. However, variable-length data is lazily checked.
