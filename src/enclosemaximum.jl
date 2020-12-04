function evalmaximum(f, intervals::Vector{Tuple{arb, arb}};
                     evaltype = :ball,
                     n = 4,
                     absmax = false)
    evals = Vector{Tuple{arb, arb}}(undef, length(intervals))
    if evaltype == :ball
        maybeabs = ifelse(absmax, abs, identity)
        for i in 1:length(intervals)
            evals[i] = getinterval(maybeabs(f(setinterval(intervals[i]...))))
        end
    elseif evaltype == :taylor
        for i in 1:length(intervals)
            evals[i] = getinterval(maximumtaylor(f, intervals[i], n, absmax = absmax))
        end
    end

    evals
end

function enclosemaximum(f,
                        a::arb,
                        b::arb;
                        # A priori computed lower bound of maximum
                        lower_bound::arb = neginf(a),
                        absmax = false,
                        evaltype = :ball,
                        n = 4,
                        # FIXME: This should be epsilon in the precision of the input parameters
                        atol = sqrt(eps()),
                        rtol = sqrt(eps()),
                        maxevals = 10^3,
                        store_trace = false,
                        show_trace = false,
                        callback = nothing,
                        extended_trace = false)
    @assert isfinite(a) && isfinite(b) && a < b

    intervals = [(a, b)]
    maxenclosure = setinterval(neginf(a), posinf(a))
    E = radius(maxenclosure)
    numevals = 0
    iteration = 0
    imprecise_input_last_iteration = false
    imprecise_input = false

    trace = MaximumTrace()
    if show_trace
        print(trace)
    end

    while (!(E < atol)
           && !(E/abs(maxenclosure) <= rtol)
           && numevals <= maxevals
           && !imprecise_input)
        iteration += 1

        # Lower and upper bound maximum on each interval
        evals = evalmaximum(f, intervals, evaltype = evaltype, n = n, absmax = absmax)
        numevals += length(intervals)

        # Find global lower and upper bound for maximum
        maxlower = neginf(a)
        maxupper = neginf(a)
        for (lower, upper) in evals
            if isfinite(lower) && !(lower < maxlower)
                maxlower = max(maxlower, lower)
            end

            if !(upper < maxupper)
                maxupper = max(maxupper, upper)
            end
        end

        if maxupper < lower_bound
            @error "given lower bound is larger than maximum"
            return parent(a)(NaN)
        end

        maxenclosure = setinterval(maxlower, maxupper)

        maxlower = max(maxlower, lower_bound)

        # Split intervals where the maximum could be located
        nextintervals = Vector{eltype(intervals)}()
        for i in 1:length(intervals)
            if !(evals[i][2] < maxlower)
                midpoint = 0.5*(intervals[i][1] + intervals[i][2])
                push!(nextintervals, (intervals[i][1], midpoint))
                push!(nextintervals, (midpoint, intervals[i][2]))
            end
        end

        # Update trace
        if !extended_trace
            dt = Dict()
        else
            dt = Dict("intervals" => intervals, "evaluations" => evals)
        end
        update!(trace,
                MaximumState(iteration, length(intervals), maxenclosure, dt),
                store_trace,
                show_trace,
                callback)

        intervals = nextintervals

        # If using Taylor series and no substantial improvement is
        # made from previous iteration this usually hints at working
        # with to low precision
        if evaltype == :taylor && isfinite(maxenclosure) && radius(maxenclosure) > 0.8*E
            if imprecise_input_last_iteration
                imprecise_input = true
            end
            imprecise_input_last_iteration = true
        else
            imprecise_input_last_iteration = false
        end

        E = radius(maxenclosure)
    end

    if imprecise_input && !(E < atol || E/abs(maxenclosure) <= rtol || !(numevals <= maxevals))
        @warn "Maximum likely needs to be computed with higher precision than $(precision(parent(a)))"
    end

    if store_trace
        return maxenclosure, trace
    end

    maxenclosure
end
