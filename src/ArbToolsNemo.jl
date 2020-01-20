module ArbToolsNemo

using Nemo
using Printf
using RecipesBase

ArbReal = arb

import Base: length, show, inv, sqrt, log, log1p, atan, asin, acos,
    exp, sin, cos, tan, sinpi, cospi, sinh, cosh, sinc, +, -, *, /, ^,
    div, sincos, isnan

import Nemo: rsqrt, tanpi, gamma, lgamma, rgamma, digamma, compose,
    derivative, integral

export arb_series

export enclosemaximum
export maximumpoly

export isolateroots

include("arb.jl")
include("arb_series.jl")

include("trace.jl")
include("enclosemaximum.jl")
include("maximumtaylor.jl")
include("maximumpoly.jl")

include("isolateroots.jl")

end # module
