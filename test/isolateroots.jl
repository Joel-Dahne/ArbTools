@testset "isolateroots" begin
    prec = 128
    RR = RealField(prec)

    problems = [(cos, 0, 10, [RR(π)/2, 3RR(π)/2, 5RR(π)/2]),
                (x -> (x - 1)*(x - 2)*(x - 3), 0, 5, [1, 2, 3]),
                (x -> exp(x) - 1, -eps(), 5, [0])]

    for (f, a, b, roots) in problems
        a = RR(a)
        b = RR(b)
        roots = RR.(roots)

        for evaltype in (:ball, :taylor)
            found, flags = isolateroots(f, a, b, evaltype = evaltype)

            @test length(found) == length(roots) && length(flags) == length(roots)
            if evaltype == :ball
                @test !any(flags)
            else
                @test all(flags)
            end
            @test all(contains(ArbTools.setinterval(found[i]...), roots[i]) for i in 1:length(roots))
        end
    end
end
