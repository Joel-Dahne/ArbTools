"""
ball: y = f(x::ArbReal)::ArbReal

taylor: f!(y::Ptr{ArbReal}, x::ArbReal, n::Integer)

"""
function evalmaximum(f, intervals::Vector{Tuple{T, T}};
                     evaltype = :ball,
                     n = 4,
                     absmax = false) where {T}
    evals = Vector{Tuple{T, T}}(undef, length(intervals))
    if evaltype == :ball
        maybeabs = ifelse(absmax, abs, identity)
        for i in 1:length(intervals)
            evals[i] = interval(maybeabs(f(setinterval(intervals[i]...))))
        end
    elseif evaltype == :taylor
        for i in 1:length(intervals)
            evals[i] = interval(maximumtaylor(f, intervals[i], n, absmax = absmax))
        end
    end

    evals
end

function enclosemaximum(f,
                        a::T,
                        b::T;
                        absmax = false,
                        evaltype = :ball,
                        n = 4,
                        atol = sqrt(eps(T)),
                        rtol = sqrt(eps(T)),
                        maxevals = 10^3,
                        store_trace = false,
                        show_trace = false,
                        extended_trace = false) where {T <: ArbReal}
    @assert isfinite(a) && isfinite(b) && a < b

    intervals = [(a, b)]
    maxenclosure = setinterval(-inf(T), inf(T))
    E = radius(maxenclosure)
    numevals = 0
    iteration = 0

    trace = MaximumTrace()
    if show_trace
        print(trace)
    end

    onehalf = ArbReal(0.5)

    # FIXME: The implementation of <= and >= in ArbNumerics are not
    # correct. If they are fixed then change < to <=.
    while !(E < atol) && !(E/abs(maxenclosure) < rtol) && numevals <= maxevals
        iteration += 1

        # Lower and upper bound maximum on each interval
        evals = evalmaximum(f, intervals, evaltype = evaltype, n = n, absmax = absmax)
        numevals += length(intervals)

        # Find global lower and upper bound for maximum
        maxlower = -inf(T)
        maxupper = -inf(T)
        for (lower, upper) in evals
            if isfinite(lower) && !(lower < maxlower)
                maxlower = max(maxlower, lower)
            end

            if !(upper < maxupper)
                maxupper = max(maxupper, upper)
            end
        end

        maxenclosure = setinterval(maxlower, maxupper)

        # Split intervals where the maximum could be located
        nextintervals = Vector{eltype(intervals)}()
        for i in 1:length(intervals)
            if !(evals[i][2] < maxlower)
                midpoint = onehalf*(intervals[i][1] + intervals[i][2])
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
                show_trace)

        intervals = nextintervals
        E = radius(maxenclosure)
    end

    #println("$n: $numevals    $maxenclosure \n\n\n\n\n\n\n\n\n\n\n\n")

    if store_trace
        return maxenclosure, trace
    end

    maxenclosure
end
