export arb_series

mutable struct arb_series
    poly::arb_poly
    length::Int
end

###############################################################################
#
#   Constructions
#
###############################################################################

function arb_series(poly::arb_poly)
    arb_series(poly, length(poly))
end

function arb_series(poly)
    arb_series(poly, length(poly))
end

function show(io::IO, x::arb_series)
    len = length(x)
    S = var(parent(x.poly))

    for i in 0:len-1
        if !iszero(x[i])
            if !isone(x[i]) || i == 0
                print(io, x[i])
                print(io, "*")
            end
            if i != 0
                print(io, "$S^$i")
            end
            print(io, " + ")
        end
    end
    print(io, "𝒪($S^$len)")
end

###############################################################################
#
#   Array and iterator interface
#
###############################################################################

function Base.getindex(x::arb_series, i::Int)
    0 <= i <= x.length || throw(BoundsError(x, i))
    return coeff(x.poly, i)
end

Base.getindex(x::arb_series, i::Number) = x[convert(Int, i)]
Base.getindex(x::arb_series, I) = [x[i] for i in I]

Base.eltype(x::arb_series) = arb
Base.iterate(x::arb_series, i = 0) = i > x.length ? nothing : (x[i], i + 1)
Base.indexed_iterate(p::arb_series, i::Int, state=1) = (x[i], i + 1)
Base.broadcastable(x::arb_series) = Ref(x)

function Base.setindex!(x::arb_series, v::arb, i::Int)
    0 <= i <= x.length || throw(BoundsError(x, i))
    setcoeff!(x.poly, i, v)
end

Base.firstindex(x::arb_series) = 0
Base.lastindex(x::arb_series) = x.length - 1

# Basic information
length(x::arb_series) = x.length

# Unary series methods
unary_methods = [# Arithmetic
                 (:inv, "arb_poly_inv_series"),
                 # Transforms
                 (:binomial_transform, "arb_poly_binomial_transform"),
                 # Powers and elementary functions
                 (:sqrt, "arb_poly_sqrt_series"),
                 (:rsqrt, "arb_poly_rsqrt_series"),
                 (:log, "arb_poly_log_series"),
                 (:log1p, "arb_poly_log1p_series"),
                 (:atan, "arb_poly_atan_series"),
                 (:asin, "arb_poly_asin_series"),
                 (:acos, "arb_poly_acos_series"),
                 (:exp, "arb_poly_exp_series"),
                 (:sin, "arb_poly_sin_series"),
                 (:cos, "arb_poly_cos_series"),
                 (:tan, "arb_poly_tan_series"),
                 (:sinpi, "arb_poly_sin_pi_series"),
                 (:cospi, "arb_poly_cos_pi_series"),
                 (:tanpi, "arb_poly_tan_pi_series"),
                 (:sinh, "arb_poly_sinh_series"),
                 (:cosh, "arb_poly_cosh_series"),
                 (:sinc, "arb_poly_sinc_series"),
                 (:sincpi, "arb_poly_sinc_pi_series"),
                 # Gamma function and factorials
                 (:gamma, "arb_poly_gamma_series"),
                 (:rgamma, "arb_poly_rgamma_series"),
                 (:lgamma, "arb_poly_lgamma_series"),
                 (:digamma, "arb_poly_digamma_series"),
                 # Zeta function
                 (:riemann_siegel_theta, "arb_poly_riemann_siegel_theta_series"),
                 (:riemann_siegel_z, "arb_poly_riemann_siegel_z_series"),
                 ]

for (s, f) in unary_methods
    @eval begin
        function ($s)(x::arb_series,
                      n::Integer = length(x))
            y = arb_series(parent(x.poly)(), n)
            ccall(($f, Nemo.libarb), Nothing,
                  (Ref{arb_poly}, Ref{arb_poly}, Int, Int),
                  y.poly, x.poly, n, precision(parent(x.poly)))
            return y
        end
    end
end

# Operators
operator_methods = [# Arithmetic
                    (:+, "arb_poly_add_series"),
                    (:-, "arb_poly_sub_series"),
                    (:*, "arb_poly_mullow"),
                    (:/, "arb_poly_div_series"),
                    # Powers and elementary functions
                    (:^, "arb_poly_pow_series"),
                    ]
for (s, f) in operator_methods
    @eval begin
        function ($s)(x::arb_series,
                      y::arb_series)
            n = max(length(x), length(y))
            z = arb_series(parent(x.poly)(), n)
            ccall(($f, Nemo.libarb), Nothing,
                  (Ref{arb_poly}, Ref{arb_poly}, Ref{arb_poly}, Int, Int),
                  z.poly, x.poly, y.poly, n, precision(parent(x.poly)))
            return z
        end
    end
end

# Binary series methods
binary_methods = [# Arithmetic
                  (:add, "arb_poly_add_series"),
                  (:sub, "arb_poly_sub_series"),
                  (:mul, "arb_poly_mullow"),
                  (:div, "arb_poly_div_series"),
                  # Composition
                  (:compose, "arb_poly_compose_series"),
                  (:revert, "arb_poly_rever_series"),
                  # Powers and elementary functions
                  (:pow, "arb_poly_pow_series"),
                  ]
for (s, f) in binary_methods
    @eval begin
        function ($s)(x::arb_series,
                      y::arb_series,
                      n::Integer = max(length(x), length(y)))
            z = arb_series(parent(x.poly)(), n)
            ccall(($f, Nemo.libarb), Nothing,
                  (Ref{arb_poly}, Ref{arb_poly}, Ref{arb_poly}, Int, Int),
                  z.poly, x.poly, y.poly, n, precision(parent(x.poly)))
            return z
        end
    end
end

# Unary length-preserving non-series methods
unary_length_preserving_methods = [(:-, "arb_poly_neg"),
                                   (:borel_transform, "arb_poly_borel_transform"),
                                   (:inv_borel_transform, "arb_poly_inv_borel_transform"),
                                   ]

for (s, f) in unary_length_preserving_methods
    @eval begin
        function ($s)(x::arb_series)
            y = arb_series(parent(x.poly)(), length(x))
            ccall(($f, Nemo.libarb), Nothing,
                  (Ref{arb_poly}, Ref{arb_poly}, Int),
                  y.poly, x.poly, precision(parent(x.poly)))
            return y
        end
    end
end

# Unary non-length-preserving methods
unary_non_length_preserving_methods = [(:derivative, "arb_poly_derivative"),
                                       (:integral, "arb_poly_integral"),
                                       ]

function derivative(x::arb_series,
                    n::Integer = max(length(x) - 1, 0))
    y = arb_series(parent(x.poly)(), n)
    ccall(("arb_poly_derivative", Nemo.libarb), Nothing,
          (Ref{arb_poly}, Ref{arb_poly}, Int),
          y.poly, x.poly, precision(parent(x.poly)))
    return y
end

function integral(x::arb_series,
                  n::Integer = length(x) + 1)
    y = arb_series(parent(x.poly)(), n)
    ccall(("arb_poly_integral", Nemo.libarb), Nothing,
          (Ref{arb_poly}, Ref{arb_poly}, Int),
          y.poly, x.poly, precision(parent(x.poly)))
    return y
end

# Special
(:^, "arb_poly_pow_ui_trunc_binexp")
(:^, "arb_poly_pow_arb_series")

(:sincos, "arb_poly_sin_cos_series")
function sincos(x::arb_series,
                n::Integer = length(x))
    y = arb_series(parent(x.poly)(), n)
    z = arb_series(parent(x.poly)(), n)
    ccall((:arb_poly_sin_cos_series, Nemo.libarb), Nothing,
          (Ref{arb_poly}, Ref{arb_poly}, Ref{arb_poly}, Int, Int),
          y.poly, z.poly, x.poly, n, precision(parent(x.poly)))
    return (y, z)
end

function atan(y::arb_series,
              x::arb_series,
              n::Integer = max(length(x), length(y)))
    xdy = mul(x, derivative(y), n - 1)
    ydx = mul(y, derivative(x), n - 1)
    x2 = mul(x, x, n - 1)
    y2 = mul(y, y, n - 1)

    z = integral((xdy - ydx)/(x2 + y2), n)
    if n > 0
        z[0] = atan(y[0], x[0])
    end
    return z
end

function zeta(s::arb_series,
              a::arb,
              n::Integer = length(s))
    y = arb_series(parent(s.poly)(), n)
    ccall((:arb_poly_zeta_series, Nemo.libarb), Nothing,
          (Ref{arb_poly}, Ref{arb_poly}, Ref{arb}, Int32, Int, Int),
          y.poly, s.poly, a, 0, n, precision(parent(s.poly)))
    return y
end
                 (:zeta, "arb_poly_zeta_series")

###############################################################################
#
#   Promotions
#
###############################################################################

promote_rule(::Type{arb_series}, ::Type{Float64}) = arb_series

promote_rule(::Type{arb_series}, ::Type{BigFloat}) = arb_series

promote_rule(::Type{arb_series}, ::Type{fmpz}) = arb_series

promote_rule(::Type{arb_series}, ::Type{fmpq}) = arb_series

promote_rule(::Type{arb_series}, ::Type{arb}) = arb_series

promote_rule(::Type{arb_series}, ::Type{fmpz_poly}) = arb_series

promote_rule(::Type{arb_series}, ::Type{fmpq_poly}) = arb_series

promote_rule(::Type{arb_series}, ::Type{T}) where {T <: Integer} = arb_series

promote_rule(::Type{arb_series}, ::Type{Rational{T}}) where T <: Union{Int, BigInt} = arb_series

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

for T in [Integer, fmpz, fmpq, Float64, BigFloat, arb, fmpz_poly, fmpq_poly, arb_poly]
    @eval begin
        +(x::arb_series, y::$T) = x + arb_series(parent(x.poly)(y))

        +(x::$T, y::arb_series) = y + x

        -(x::arb_series, y::$T) = x - arb_series(parent(x.poly)(y))

        -(x::$T, y::arb_series) = arb_series(parent(y.poly)(x)) - y

        *(x::arb_series, y::$T) = x * arb_series(parent(x.poly)(y))

        *(x::$T, y::arb_series) = y * x

        /(x::arb_series, y::$T) = x / arb_series(parent(x.poly)(y))

        /(x::$T, y::arb_series) = arb_series(parent(y.poly)(x)) / y

        ^(x::arb_series, y::$T) = x^arb_series(parent(x.poly)(y))

        ^(x::$T, y::arb_series) = arb_series(parent(y.poly)(x))^y
    end
end

+(x::arb_series, y::Rational{T}) where T <: Union{Int, BigInt} = x + arb_series(parent(x.poly)(y))

+(x::Rational{T}, y::arb_series) where T <: Union{Int, BigInt} = y + x

-(x::arb_series, y::Rational{T}) where T <: Union{Int, BigInt} = x - arb_series(parent(x.poly)(y))

-(x::Rational{T}, y::arb_series) where T <: Union{Int, BigInt} = arb_series(parent(y.poly)(x)) - y

*(x::arb_series, y::Rational{T}) where T <: Union{Int, BigInt} = x * arb_series(parent(x.poly)(y))

*(x::Rational{T}, y::arb_series) where T <: Union{Int, BigInt} = y * x

/(x::arb_series, y::Rational{T}) where T <: Union{Int, BigInt} = x / arb_series(parent(x.poly)(y))

/(x::Rational{T}, y::arb_series) where T <: Union{Int, BigInt} = arb_series(parent(y.poly)(x)) / y

^(x::arb_series, y::Rational{T}) where T <: Union{Int, BigInt} = x^arb_series(parent(x.poly)(y))

^(x::Rational{T}, y::arb_series) where T <: Union{Int, BigInt} = arb_series(parent(y.poly)(x))^y
