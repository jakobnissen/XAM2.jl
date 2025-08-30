module SAMFlag
    export Flag, PAIRED, PROPER_PAIR, UNMAP, MUNMAP, REVERSE, MREVERSE, READ1, READ2, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY

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

struct SAMFlags <: AbstractSet{Flag}
    x::UInt16

    global new_sam_flags(x::UInt16) = new(x)
end

SAMFlags(x::UInt16) = new_sam_flags(x & 0x0fff)

Base.Integer(x::SAMFlags) = x.x

function Base.iterate(x::SAMFlags, state::UInt16 = x.x)
    return if iszero(state)
        nothing
    else
        (reinterpret(Flag, trailing_zeros(state) % UInt8), state & (state - 0x0001))
    end
end

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
