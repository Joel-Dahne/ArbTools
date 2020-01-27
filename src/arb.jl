function Base.eps(x::arb)
    parent(x)(2)^(-(prec(parent(x)) - 1))
end

function Base.eps(RR::ArbField)
    RR(2)^(-(prec(RR) - 1))
end

function isnan(x::arb)
    x_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), x)
    0 != ccall(("arf_is_nan", Nemo.libarb), Cint, (Ref{Nemo.arf_struct},), x_mid)
end

function inf(x::arb)
    y = parent(x)(0)
    ccall(("arb_pos_inf", Nemo.libarb), Cvoid, (Ref{arb},), y)
    y
end

function posinf(x::arb)
    y = parent(x)(0)
    ccall(("arb_pos_inf", Nemo.libarb), Cvoid, (Ref{arb},), y)
    y
end

function neginf(x::arb)
    y = parent(x)(0)
    ccall(("arb_neg_inf", Nemo.libarb), Cvoid, (Ref{arb},), y)
    y
end

"""
    max(x::arb, y::arb)
> Return a ball containing the maximum of x and y.
"""
function Base.max(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_max, Nemo.libarb), Nothing,
          (Ref{arb}, Ref{arb}, Ref{arb}, Int), z, x, y, parent(x).prec)
    return z
end

"""
    min(x::arb, y::arb)
> Return a ball containing the minimum of x and y.
"""
function Base.min(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_min, Nemo.libarb), Nothing,
          (Ref{arb}, Ref{arb}, Ref{arb}, Int), z, x, y, parent(x).prec)
    return z
end

"""
    setinterval(x::arb, y::arb)
> Return a ball containing the interval [x,y].
"""
function setinterval(x::arb, y::arb)
    z = parent(x)()
    x_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), x)
    y_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), y)
    ccall((:arb_set_interval_arf, Nemo.libarb), Cvoid,
          (Ref{arb}, Ptr{Nemo.arf_struct}, Ptr{Nemo.arf_struct}, Int),
          z, x_mid, y_mid, x.parent.prec)
    return z
end

"""
    getinterval(x::arb)
    getinterval(::Type{arb}, x::arb)
> Return an interval [a,b] containing the ball x.
"""
function getinterval(x::arb)
    getinterval(arb, x)
end

function getinterval(::Type{arb}, x::arb)
    a, b = x.parent(), x.parent()
    a_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), a)
    b_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), b)
    ccall((:arb_get_interval_arf, Nemo.libarb), Cvoid,
          (Ptr{Nemo.arf_struct}, Ptr{Nemo.arf_struct}, Ref{arb}, Clong),
          a_mid, b_mid, x, x.parent.prec)

    (a, b)
end

"""
    getinterval(::Type{BigFloat}, x::arb)
> Return an interval [a,b] containing the ball x.
"""
function getinterval(::Type{BigFloat}, x::arb)
    a, b = BigFloat(), BigFloat()
    ccall((:arb_get_interval_mpfr, Nemo.libarb), Cvoid,
          (Ref{BigFloat}, Ref{BigFloat}, Ref{arb}),
          a, b, x)

    (a, b)
end

"""
    convert(::Type{BigFloat}, x::arb)
> Return the midpoint of x as a BigFloat rounded down to the current
  precision of BigFloat.
"""
function BigFloat(x::arb)
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), x)
        # 4 == round to nearest
        m = BigFloat()
        ccall((:arf_get_mpfr, Nemo.libarb), Float64,
              (Ref{BigFloat}, Ptr{Nemo.arf_struct}, Base.MPFR.MPFRRoundingMode),
              m, t, Base.MPFR.MPFRRoundNearest)
    end
    return m
end

"""
    convert(::Type{BigFloat}, x::arb)
> Return the midpoint of x as a BigFloat rounded down to the current
  precision of BigFloat.
"""
function Base.convert(::Type{BigFloat}, x::arb)
    return BigFloat(x)
end

"""
    rel_accuracy_bits(x::arb)
> Compute the relatively accuracy of the ball `x` in bits.
"""
function rel_accuracy_bits(x::arb)
    ccall(("arb_rel_accuracy_bits", Nemo.libarb), Int,
          (Ref{arb},), x)
end

"""
    atan(x::arb, y::arb)
> Return atan(x, y) = arg(x + yi).
"""
function atan(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_atan2, Nemo.libarb), Nothing,
          (Ref{arb}, Ref{arb}, Ref{arb}, Int), z, x, y, parent(x).prec)
    return z
end
