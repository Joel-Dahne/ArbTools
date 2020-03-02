# ArbTools.jl

This package implements some tools for working with rigorous numerics.
It relies on the [Arb C Library](http://arblib.org/index.html) for all
computations and uses [Nemo](https://github.com/wbhart/Nemo.jl) as the
basic interface between Arb and Julia.

The package adds some extra functionality to the `arb` type from Nemo,
it defines a new type `arb_series` for working with Taylor series of
functions and two methods, `enclosemaximum` for rigorously bounding
the maximum value of a univariate function and `isolateroots` to
rigorously isolate all roots of a univarate function.

This package is mainly intended to implement functionality that I need
in my research but that is missing in Nemo. It should not be
considered a stable package and will change as my needs change.

## Installation
The package is not in the general Julia repository but can be
installed through the package manager with
``` julia
pkg> add https://github.com/Joel-Dahne/ArbTools.jl
```
