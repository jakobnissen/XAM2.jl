"""
    module SAMFlag

This module works as a namespace for instances of the enum SAM flags.
You can access a given SAM flag with e.g. `SAMFlag.SUPPLEMENTARY`.
"""
module SAMFlag
    export Flag, PAIRED, PROPER_PAIR, UNMAP, MUNMAP, REVERSE, MREVERSE, READ1, READ2, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY

    """
        Flag

    A 1-byte primitive type representing a SAM flag. Its value is guaranteed to
    be equal to the number of trailing byte in the SAM binary representation,
    e.g. in SAM, `UNMAP` is 4, so
    `reinterpret(UInt8, SAMFlag.UNMAP) == trailing_zeros(4)`.

    Access `Flag`s through the `SAMFlag` module namespace.

    The following `Flag`s currenlty exists. More may be added as the SAM
    spefication evolves:
    `PAIRED`, `PROPER_PAIR`, `UNMAP`, `MUNMAP`, `REVERSE`, `MREVERSE`,
    `READ1`, `READ2`, `SECONDARY`, `QCFAIL`, `DUP`, `SUPPLEMENTARY`
    """
    @enum Flag::UInt8 begin
        PAIRED
        PROPER_PAIR
        UNMAP
        MUNMAP
        REVERSE
        MREVERSE
        READ1
        READ2
        SECONDARY
        QCFAIL
        DUP
        SUPPLEMENTARY
    end
end

using .SAMFlag

"""
    SAMFlags <: AbstractSet{Flag}

A compact, immutable, efficient set of `Flag`. `Integer(x::SAMFlags)` will return a
`UInt16` equal to their representation in the SAM format, e.g. the set
of flags `PAIRED, UNMAP` is encoded as `5` is SAM, and so
`Integer(SAMFlags([SAMFlag.PAIRED, SAMFlag.UNMAP])) === UInt16(5)`.

Construct from a `UInt16` with its integer representation, or an iterable
of `Flag`.

# Examples
```jldoctest
julia> flags = SAMFlags(0x0005)
SAMFlags with 2 elements:
  XAM2.SAMFlag.PAIRED
  XAM2.SAMFlag.UNMAP

julia> Integer(flags)
0x0005

julia> flags === SAMFlags([SAMFlag.PAIRED, SAMFlag.UNMAP])
true

julia> SAMFlag.PAIRED in flags
true

julia> length(intersect(flags, SAMFlags([SAMFlag.UNMAP])))
1
```
"""
struct SAMFlags <: AbstractSet{Flag}
    x::UInt16

    global new_sam_flags(x::UInt16) = new(x)
end

SAMFlags(x::UInt16) = new_sam_flags(x & 0x0fff)

function SAMFlags(itr)
    flags = SAMFlags(0x0000)
    for i in itr
        flag = convert(Flag, i)::Flag
        flags = push(flags, flag)
    end
    return flags
end

Base.Integer(x::SAMFlags) = x.x

function Base.iterate(x::SAMFlags, state::UInt16 = x.x)
    return if iszero(state)
        nothing
    else
        (reinterpret(Flag, trailing_zeros(state) % UInt8), state & (state - 0x0001))
    end
end

"""
    push(flags::SAMFlags, flag::Flag)

Return a new `SAMFlags` containing `flag` and all flags in `flags`.

# Example
```
julia> flags = SAMFlags([SAMFlag.DUP, SAMFlag.READ1]);

julia> flags2 = push(flags, SAMFlag.PAIRED)
SAMFlags with 3 elements:
  XAM2.SAMFlag.PAIRED
  XAM2.SAMFlag.READ1
  XAM2.SAMFlag.DUP

julia> length(flags) # old one unchanged
2
```
"""
push(x::SAMFlags, y::Flag) = new_sam_flags(x.x | (UInt16(1) << (reinterpret(UInt8, y) & 15)))
Base.in(flag::Flag, x::SAMFlags) = isodd(x.x >>> (reinterpret(UInt8, flag) & UInt8(15)))
Base.length(x::SAMFlags) = count_ones(x.x)
Base.isempty(x::SAMFlags) = iszero(x.x)
Base.hash(x::SAMFlags, h::UInt) = hash(x.x, h ⊻ (0x94e973b118a1b3ae % UInt))
Base.union(x::SAMFlags, y::SAMFlags) = new_sam_flags(x.x | y.x)
Base.intersect(x::SAMFlags, y::SAMFlags) = new_sam_flags(x.x & y.x)
Base.setdiff(x::SAMFlags, y::SAMFlags) = new_sam_flags(x.x & ~y.x)
Base.issubset(x::SAMFlags, y::SAMFlags) = isempty(setdiff(x, y))
Base.isdisjoint(x::SAMFlags, y::SAMFlags) = isempty(x ∩ y)
