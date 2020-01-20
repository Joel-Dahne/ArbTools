function isuniqueroot(f,
                      a::arb,
                      b::arb;
                      evaltype = :taylor)
    if evaltype == :ball
        @warn ("Cannot determine if a root is unique with only "
               *"evaluation of f, it requires an evaluation of the Taylor "
               *"expansion")

        return false, true
    elseif evaltype == :taylor
        x = setinterval(a, b)
        PP = ArbPolyRing(parent(x), :x)

        # Evaluate f' on the interval and check that it's non-zero
        dy = f(arb_series(PP([x, parent(x)(1)])))

        if contains_zero(dy[1])
            return false, true
        end

        # Evaluate f at the endpoints and check that they have different signs
        ya = f(a)
        yb = f(b)

        if ispositive(ya) && isnegative(yb) || isnegative(ya) && ispositive(yb)
            return true, true
        end

        if ispositive(ya) && ispositive(yb) || isnegative(ya) && isnegative(yb)
            return false, false
        end

        return false, true
    end

    return false, true
end

function refine_root(poly::arb_poly,
                     (a, b)::Tuple{arb, arb})
    root = setinterval(a, b)
    # Perform a couple of simple newton iterations to gain at least 20
    # bits of relative accuracy
    iterations = 1
    while iterations < 5 && rel_accuracy_bits(root) < 20
        iterations += 1
        y = evaluate(poly, midpoint(root))
        dy = evaluate(derivative(poly), root)

        root = setintersection(midpoint(root) - y/dy, root)

        if isnan(root)
            return (a, b)
        end
    end

    # Refine the root using method in arb
    convergence_interval = parent(root)(root)
    convergence_factor_ball = parent(root)(0)

    GC.@preserve convergence_factor_ball begin
        convergence_factor = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct},
                                   (Ref{arb}, ), convergence_factor_ball)

        ccall(("_arb_poly_newton_convergence_factor", Nemo.libarb), Cvoid,
              (Ptr{Nemo.arf_struct}, Ptr{arb}, Int, Ref{arb}, Int),
              convergence_factor, poly.coeffs, length(poly), convergence_interval, parent(root).prec)

        # FIXME: The extra bits used in the computation should be
        # looked over. I just chose it randomly.
        ccall(("_arb_poly_newton_refine_root", Nemo.libarb), Cvoid,
              (Ref{arb}, Ptr{arb}, Int, Ref{arb}, Ref{arb}, Ref{Nemo.arf_struct}, Int, Int),
              root, poly.coeffs, length(poly), root, convergence_interval,
              convergence_factor, div(parent(root).prec, 4), parent(root).prec)
    end

    if a < root
        return getinterval(root)
    else
        return (a, b)
    end
end

function isolateroots(f,
                      a::arb,
                      b::arb;
                      evaltype = :ball,
                      refine = false,
                      atol = sqrt(eps()),
                      rtol = sqrt(eps()),
                      maxevals = 10^3,
                      maxfound = typemax(Int),
                      store_trace = false,
                      show_trace = false,
                      extended_trace = false)
    @assert isfinite(a) && isfinite(b) && a < b

    intervals = [(a, b)]
    numfound = 0
    numevals = 0
    iteration = 0

    found = Vector{eltype(intervals)}()
    flags = Vector{Bool}()

    while !isempty(intervals) && numfound < maxfound && numevals < maxevals
        iteration += 1
        nextintervals = Vector{eltype(intervals)}()

        for (lower, upper) in intervals
            numevals += 1

            mayberoot = contains_zero(f(setinterval(lower, upper)))
            unique = false

            if mayberoot
                if evaltype == :taylor
                    # Try to prove uniqueness
                    unique, mayberoot = isuniqueroot(f, lower, upper, evaltype = evaltype)
                end

                if unique && refine
                    # Try to refine the root
                    (lower, upper) = refine_root(f, lower, upper,
                                                 evaltype = evaltype,
                                                 atol = atol,
                                                 rtol = rtol)
                end
            end

            if mayberoot
                x = setinterval(lower, upper)

                if (radius(x) < atol) || (radius(x)/x < rtol)
                    push!(found, (lower, upper))
                    push!(flags, unique)
                else
                    midpoint = 0.5*(lower + upper)
                    push!(nextintervals, (lower, midpoint))
                    push!(nextintervals, (midpoint, upper))
                end
            end
        end

        if show_trace
            println("$iteration      $(length(nextintervals))")
        end

        intervals = nextintervals
    end

    found = [found; intervals]
    flags = [flags; zeros(Bool, length(intervals))]

    p = sortperm(found, by=x -> BigFloat(getindex(x, 1)))

    found[p], flags[p]
end

function isolateroots(poly::arb_poly,
                      a::arb,
                      b::arb;
                      atol = sqrt(eps()),
                      rtol = sqrt(eps()),
                      maxevals = 10^3,
                      maxfound = typemax(Int),
                      store_trace = false,
                      show_trace = false,
                      extended_trace = false)
    @assert isfinite(a) && isfinite(b) && a < b

    intervals = [(a, b)]
    numfound = 0
    numevals = 0
    iteration = 0

    found = Vector{eltype(intervals)}()
    flags = Vector{Bool}()

    while !isempty(intervals) && numfound < maxfound && numevals < maxevals
        iteration += 1
        nextintervals = Vector{eltype(intervals)}()

        for (lower, upper) in intervals
            numevals += 1

            x = setinterval(lower, upper)
            y = evaluate(poly, x)

            if contains_zero(y)
                dy = evaluate(derivative(poly), x)

                if contains_zero(dy)
                    ylower = evaluate(poly, lower)
                    yupper = evaluate(poly, upper)

                    if ylower > 0 && yupper < 0 || ylower < 0 && yupper > 0
                        (lower, upper) = refine_root(poly, (lower, upper))
                        push!(found, (lower, upper))
                        push!(flags, true)
                    elseif !(ylower > 0 && yupper > 0 || ylower < 0 && yupper < 0) &&
                        if (radius(x) < atol) || (radius(x)/x < rtol)
                            push!(found, (lower, upper))
                            push!(flags, false)
                        else
                            midpoint = 0.5*(lower + upper)
                            push!(nextintervals, (lower, midpoint))
                            push!(nextintervals, (midpoint, upper))
                        end
                    end
                elseif (radius(x) < atol) || (radius(x)/x < rtol)
                    push!(found, (lower, upper))
                    push!(flags, false)
                else
                    midpoint = 0.5*(lower + upper)
                    push!(nextintervals, (lower, midpoint))
                    push!(nextintervals, (midpoint, upper))
                end
            end
        end

        if show_trace
            println("$iteration      $(length(nextintervals))")
        end

        intervals = nextintervals
    end


    found = [found; intervals]
    flags = [flags; zeros(Bool, length(intervals))]

    p = sortperm(found, by=x -> BigFloat(getindex(x, 1)))

    found[p], flags[p]
end
