function maximumtaylor(f!, xinterval, n; absmax = false)
    x = setinterval(xinterval...)
    maybeabs = ifelse(absmax, abs, identity)

    poly = _arb_vec_init(n)

    # Compute the Taylor polynomial at the midpoint of x
    f!(poly, midpoint(x), n)

    # The maximum is lower bounded by the value at the midpoint
    #lower = lowerbound(maybeabs(unsafe_load_ArbRealPtr(poly, 1)))

    # Enclose the Taylor polynomial evaluated at x - midpoint(x)
    #res = maybeabs(_arb_poly_evaluate(poly, n, x - midpoint(x)))
    res = maximumpoly(poly, n, xinterval .- midpoint(x), absmax = absmax)
    #@show res

    # Compute the rest term of the Taylor expansion
    f!(poly, x, n)

    restterm = pow((x - midpoint(x)), n - 1)*unsafe_load_ArbRealPtr(poly, n)
    #@show radius(restterm)
    res += restterm
    #upper = upperbound(res)

    _arb_vec_clear(poly, n)

    #return setinterval(lower, upper)
    return res
end
