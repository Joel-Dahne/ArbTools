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
    x_lower = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), lbound(x))
    y_upper = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), ubound(y))
    ccall((:arb_set_interval_arf, Nemo.libarb), Cvoid,
          (Ref{arb}, Ptr{Nemo.arf_struct}, Ptr{Nemo.arf_struct}, Int),
          z, x_lower, y_upper, x.parent.prec)
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
    rel_accuracy_bits(x)
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

function arb_dump(x::arb)
    cstr = ccall((:arb_dump_str, Nemo.libarb), Ptr{UInt8}, (Ref{arb},),
                 x)

    unsafe_string(cstr)
end

function arb_load_dump(str::String, r::ArbField)
    x = r()
    err = ccall((:arb_load_str, Nemo.libarb), Int32, (Ref{arb}, Ptr{UInt8}),
                x, str)
    err == 0 || Throw(error("Invalid string $str"))
    x
end

function format_arb(x::arb, digits::Int)
    cstr = ccall((:arb_get_str, :libarb), Ptr{UInt8}, (Ref{arb}, Int, UInt),
                 x, digits, UInt(0))
    str = unsafe_string(cstr)
    ccall((:flint_free, Nemo.libflint), Nothing, (Ptr{UInt8},), cstr)
    return str
end

function add_error!(x::arb, error::arb)
    ccall((:arb_add_error, Nemo.libarb), Nothing, (Ref{arb}, Ref{arb}), x, error)
    return x
end

function abs_ubound(x::arb)
    res = parent(x)()
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), res)
        ccall((:arb_get_abs_ubound_arf, Nemo.libarb), Nothing, (Ptr{Nemo.arf_struct}, Ref{arb}, Int),
              t, x, parent(x).prec)
    end
    return res
end

function abs_lbound(x::arb)
    res = parent(x)()
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), res)
        ccall((:arb_get_abs_lbound_arf, Nemo.libarb), Nothing, (Ptr{Nemo.arf_struct}, Ref{arb}, Int),
              t, x, parent(x).prec)
    end
    return res
end

function ubound(x::arb)
    res = parent(x)()
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), res)
        ccall((:arb_get_ubound_arf, Nemo.libarb), Nothing, (Ptr{Nemo.arf_struct}, Ref{arb}, Int),
              t, x, parent(x).prec)
    end
    return res
end

function lbound(x::arb)
    res = parent(x)()
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), res)
        ccall((:arb_get_lbound_arf, Nemo.libarb), Nothing, (Ptr{Nemo.arf_struct}, Ref{arb}, Int),
              t, x, parent(x).prec)
    end
    return res
end
