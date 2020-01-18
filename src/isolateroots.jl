"""
ball: y = f(x::ArbReal)::ArbReal

taylor: f!(y::Ptr{ArbReal}, x::ArbReal, n::Integer)

"""
function mayberoot(f, lower, upper;
                   evaltype = :ball)
    if evaltype == :ball
        return contains(f(setinterval(lower, upper)), zero(lower))
    elseif evaltype == :taylor
        y = _arb_vec_init(1)
        f(y, setinterval(lower, upper), 1)
        maybe = contains(unsafe_load_ArbRealPtr(y, 1), zero(lower))
        _arb_vec_clear(y, 1)
        return maybe
    end
end

function isuniqueroot(f, lower, upper;
                      evaltype = :ball)
    if evaltype == :ball
        error("Cannot determine if a root is unique with only
            evaluation of f, it requires an evaluation of the Taylor
            expansion")
    elseif evaltype == :taylor
        x = setinterval(lower, upper)
        # Evaluate f at the midpoint
        y = _arb_vec_init(1)
        f(y, midpoint(x), 1)
        fmid = unsafe_load_ArbRealPtr(y, 1)
        _arb_vec_clear(y, 1)

        # Evaluate f' on the interval
        y = _arb_vec_init(2)
        f(y, x, 2)
        df = unsafe_load_ArbRealPtr(y, 2)
        _arb_vec_clear(y, 2)

        enclosure = midpoint(x) - fmid/df

        if isnan(enclosure)
            return (lower, upper), false, true
        elseif contains(x, enclosure)
            return interval(enclosure), true, true
        elseif enclosure < x || enclosure > x
            return interval(enclosure), false, false
        else
            return interval(setintersection(x, enclosure)), false, true
        end
    end

    (lower, upper), false, true
end

function isolateroots(f,
                      a::T,
                      b::T;
                      evaltype = :ball,
                      atol = sqrt(eps(T)),
                      rtol = sqrt(eps(T)),
                      maxevals = 10^3,
                      maxfound = typemax(Int),
                      store_trace = false,
                      show_trace = false,
                      extended_trace = false) where {T <: ArbReal}
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
            if mayberoot(f, lower, upper, evaltype = evaltype)
                x = setinterval(lower, upper)

                if evaltype == :ball
                    # For evaltype = :ball we cannot hope to prove
                    # that the there is a unique root
                    if (radius(x) < atol) || (radius(x)/x < rtol)
                        push!(found, (lower, upper))
                        push!(flags, false)
                    else
                        midpoint = 0.5*(lower + upper)
                        push!(nextintervals, (lower, midpoint))
                        push!(nextintervals, (midpoint, upper))
                    end
                elseif evaltype == :taylor
                    numevals += 1
                    (lower, upper), unique, maybe = isuniqueroot(f, lower, upper, evaltype = evaltype)

                    if unique
                        push!(found, (lower, upper))
                        push!(flags, true)
                    elseif maybe && ((radius(x) < atol) || (radius(x)/x < rtol))
                        push!(found, (lower, upper))
                        push!(flags, false)
                    elseif maybe
                        midpoint = 0.5*(lower + upper)
                        push!(nextintervals, (lower, midpoint))
                        push!(nextintervals, (midpoint, upper))
                    end
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

    p = sortperm(found, by=x -> getindex(x, 1))

    found[p], flags[p]
end
