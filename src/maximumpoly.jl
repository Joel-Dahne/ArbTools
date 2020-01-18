function maximumpoly(poly::Ptr{ArbReal}, n::Int,
                     xinterval::Tuple{ArbReal, ArbReal};
                     absmax::Bool = false)
    lower, upper = xinterval
    x = setinterval(xinterval...)
    maybeabs = ifelse(absmax, abs, identity)

    # Enclose the zeros of the derivative of the polynomial
    deriv = _arb_vec_init(n)

    _arb_poly_derivative(deriv, poly, n, workingprecision(x))

    f! = (y, x, k) -> begin
        if k == 1
            unsafe_store_ArbRealPtr!(y, _arb_poly_evaluate(deriv, n, x), 1)
        elseif k == 2
            (a, b) = ArbTools._arb_poly_evaluate2(deriv, n, x)
            unsafe_store_ArbRealPtr!(y, a, 1)
            unsafe_store_ArbRealPtr!(y, b, 2)
        end
        nothing
    end

    roots = isolateroots(f!, lower, upper, evaltype = :taylor)[1]

    _arb_vec_clear(deriv, n)

    # Evaluate the polynomial on the zeros and the endpoints and take
    # the maximum
    m = max(maybeabs(_arb_poly_evaluate(poly, n, lower)),
            maybeabs(_arb_poly_evaluate(poly, n, upper)))

    for root in roots
        m = max(m,
                maybeabs(_arb_poly_evaluate(poly, n, setinterval(root...))))
    end

    m
end
