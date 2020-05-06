function maximumtaylor(f, (a, b), n; absmax = false)
    x = setinterval(a, b)
    maybeabs = ifelse(absmax, abs, identity)
    PP = ArbPolyRing(parent(x), :x)

    # Compute the rest term of the Taylor expansion
    y = f(arb_series(PP([x, parent(x)(1)]), n + 1))

    restterm = (x - midpoint(x))^n*y[n]

    # If the rest term is not finite the result always be not finite
    if !isfinite(restterm)
        return restterm
    end

    # Compute the Taylor polynomial at the midpoint of x
    y = f(arb_series(PP([midpoint(x), parent(x)(1)]), n))

    # Enclose the Taylor polynomial evaluated at (a, b) - midpoint((a, b))
    mid = 0.5*(a + b)
    res = maximumpoly(y.poly, (a - mid, b - mid), absmax = absmax)
    p = res

    res += restterm

    return res
end
