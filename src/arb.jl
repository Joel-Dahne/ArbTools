"""
This file contains functions related to Arb which should probably be
in the ArbNumerics package but currently is not.
"""

"""
    contains(x, y)

    Returns true iff the given number (or ball) y is contained in the
    interval represented by x.

    If x contains NaN, this function always returns nonzero (as it could
    represent anything, and in particular could represent all the points
    included in y). If y contains NaN and x does not, it always returns
    zero.
"""
function contains(x::ArbReal, y::ArbReal)
    0 != ccall(ArbNumerics.@libarb(arb_contains), Cint, (Ref{ArbReal}, Ref{ArbReal}), x, y)
end

"""
    overlaps(x, y)

    Returns nonzero iff x and y have some point in common. If either x
    or y contains NaN, this function always returns nonzero (as a NaN
    could be anything, it could in particular contain any number that
    is included in the other operand). """

#        int arb_overlaps(const arb_t x, const arb_t y)
"""
    setunion(x, y)

    Returns a ball containing both x and y.
"""
function setunion(x::ArbReal{P}, y::ArbReal{P}) where {P}
    z = ArbReal{P}()
    ccall(ArbNumerics.@libarb(arb_union), Cvoid,
          (Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}), z, x, y)
    return z
end

"""
    setintersection(x, y)

    If x and y overlap according to arb_overlaps(), then returns a
    ball containing the intersection of x and y. Otherwise NaN is
    returned.
"""
function setintersection(x::ArbReal{P}, y::ArbReal{P}) where {P}
    z = ArbReal{P}()
    res = ccall(ArbNumerics.@libarb(arb_intersection), Cint,
                (Ref{ArbReal}, Ref{ArbReal}, Ref{ArbReal}), z, x, y)
    if res != 0
        return z
    else
        return ArbReal{P}(NaN)
    end
end

"""
    _arb_vec_init(n)

    Returns a pointer to an array of n initialized arb_struct entries.
"""
function _arb_vec_init(n::Int)
    ccall((:_arb_vec_init, :libarb), Ptr{ArbReal}, (Int,), n)
end

"""
    _arb_vec_clear(v, n)

    Clears an array of n initialized arb_struct entries.
"""
function _arb_vec_clear(v::Ptr{ArbReal}, n::Int)
    ccall((:_arb_vec_clear, :libarb), Cvoid, (Ptr{ArbReal}, Int), v, n)
end

function unsafe_load_ArbRealPtr(ptr::Ptr{ArbReal}, i::Int)
    res = ArbReal(0)
    ccall((:arb_set, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal}),
          res, ptr + (i - 1)*sizeof(ArbReal))

    return res
end

function unsafe_store_ArbRealPtr!(ptr::Ptr{ArbReal}, value::ArbReal, i::Int)
    ccall((:arb_set, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal}),
          ptr + (i - 1)*sizeof(ArbReal), value)
end

function _arb_poly_evaluate(f::Ptr{ArbReal}, len::Integer, x::ArbReal)
    res = zero(x)
    ccall((:_arb_poly_evaluate, :libarb), Cvoid, (Ref{ArbReal}, Ptr{ArbReal}, Clong,
                                                   Ref{ArbReal}, Clong),
          res, f, len, x, workingprecision(x))

    return res
end

function _arb_poly_evaluate2(f::Ptr{ArbReal}, len::Integer, x::ArbReal)
    y = zero(x)
    z = zero(x)
    ccall((:_arb_poly_evaluate2, :libarb), Cvoid, (Ref{ArbReal}, Ref{ArbReal},
                                                   Ptr{ArbReal}, Clong,
                                                   Ref{ArbReal}, Clong),
          y, z, f, len, x, workingprecision(x))

    return (y, z)
end

function _arb_poly_derivative(res::Ptr{ArbReal}, poly::Ptr{ArbReal}, len::Int, prec::Int)
    ccall((:_arb_poly_derivative, :libarb), Cvoid, (Ptr{ArbReal}, Ptr{ArbReal},
                                                    Clong, Clong),
          res, poly, len, prec)
end
