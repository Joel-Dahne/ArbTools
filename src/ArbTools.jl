module ArbTools

using Nemo
using Printf
using RecipesBase
using Colors

ArbReal = arb

import Base: length, show, inv, sqrt, log, log1p, atan, asin, acos,
    exp, sin, cos, tan, sinpi, cospi, sinh, cosh, sinc, +, -, *, /, ^,
    div, sincos, isnan

import Nemo: rsqrt, tanpi, gamma, lgamma, rgamma, digamma, compose,
    derivative, integral, zeta

export setinterval, getinterval

export arb_series

export enclosemaximum
export maximumpoly

export isolateroots

include("arb.jl")
include("acb.jl")
include("arb_series.jl")
include("acb_series.jl")

include("trace.jl")
include("enclosemaximum.jl")
include("maximumtaylor.jl")
include("maximumpoly.jl")

include("isolateroots.jl")

include("integrate.jl")

end # module
