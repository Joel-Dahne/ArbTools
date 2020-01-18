struct MaximumState{P}
    iteration::Int
    numintervals::Int
    maxenclosure::ArbReal{P}
    metadata::Dict
end

function MaximumState(i, numintervals, maxenclosure)
    MaximumState(Int(i), numintervals, maxenclosure, Dict())
end

struct MaximumTrace
    states::Vector{MaximumState}
end

MaximumTrace() = MaximumTrace(Array{MaximumState}(undef, 0))

function Base.show(io::IO, s::MaximumState)
    @printf io "%6d   %11d    %s\n" s.iteration s.numintervals string(s.maxenclosure)
    return
end

function Base.show(io::IO, t::MaximumTrace)
    @printf io "Iter     Numintervals   Enclosure \n"
    @printf io "------   ------------   --------------\n"
    for state in t.states
        show(io, state)
    end
    return
end

Base.push!(t::MaximumTrace, s::MaximumState) = push!(t.states, s)

function update!(tr::MaximumTrace,
                 st::MaximumState,
                 store_trace::Bool,
                 show_trace::Bool)
    if store_trace
        push!(tr, st)
    end
    if show_trace
        show(st)
    end
    return
end

@recipe function f(st::MaximumState)
    xs = [Float64(interval[j]) for j in 1:2, interval in st.metadata["intervals"]]
    ys = [[Float64(eval[j]) for eval in st.metadata["evaluations"]] for j in 1:2]
    seriescolor --> :blue
    label := ""
    fillrange := ys[2]'
    (xs, ys[1]')
end

@userplot TracePlot

@recipe function f(h::TracePlot; cname = "oranges")
    if typeof(h.args[1]) == MaximumTrace
        states = h.args[1].states
    elseif typeof(h.args[1]) == Vector{MaximumState}
        states = h.args[1]
    else
        error("Trace plots should be given a trace or vector of states.  Got: $(typeof(h.args[1]))")
    end

    N = length(states)

    coloroffset = 1
    cm = colormap(cname, coloroffset + N)

    for i in 1:N
        @series begin
            seriescolor := cm[i + coloroffset]
            states[i]
        end
    end
end
