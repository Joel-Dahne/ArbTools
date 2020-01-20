function maximumpoly(poly::arb_poly,
                     (lower, upper)::Tuple{arb, arb};
                     absmax::Bool = false)
    x = setinterval(lower, upper)
    maybeabs = ifelse(absmax, abs, identity)

    # Enclose the zeros of the derivative of the polynomial
    deriv = derivative(poly)

    roots = isolateroots(deriv, lower, upper)[1]

    # Evaluate the polynomial on the zeros and the endpoints and take
    # the maximum
    m = max(maybeabs(evaluate(poly, lower)),
            maybeabs(evaluate(poly, upper)))

    # FIXME: If lower and upper are not tight then it can happen that
    # we find a root which is contained in x but not in the true
    # interval. If the maximum occurs at this root we would get wrong
    # results.
    for root in roots
        m = max(m,
                maybeabs(evaluate(poly, setinterval(root...))))
    end

    m
end
