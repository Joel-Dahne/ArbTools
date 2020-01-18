module ArbToolsNemo

export enclosemaximum
export maximumpoly

export isolateroots

using Nemo
using ArbNumerics
using Printf
using RecipesBase

include("arb.jl")

include("trace.jl")
include("enclosemaximum.jl")
include("maximumtaylor.jl")
include("maximumpoly.jl")

include("isolateroots.jl")

end # module
