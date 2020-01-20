@testset "continuous" begin
    prec = 128
    RR = RealField(prec)
    atol = sqrt(eps())
    rtol = sqrt(eps())

    problems = [(sin, 0, 1, sin(RR(1))),
                (sinpi, 0, 1, 1),
                (x -> 2x + 1, 0, 1, 3),
                (x -> -(x + 1e-1)^3 , 0, 1, 0)]

    for (f, a, b, value) in problems
        a = RR(a)
        b = RR(b)
        value = RR(value)

        for evaltype in (:ball, :taylor)
            enclosure = enclosemaximum(f, a, b, atol = atol, rtol = rtol, evaltype = evaltype)

            @test contains(enclosure, value)
            @test radius(value) < atol || radius(value)/abs(value) < rtol
        end
    end
end

@testset "unbounded" begin
    prec = 128
    RR = RealField(prec)

    a = RR(0)
    b = RR(1)

    @test isnan(enclosemaximum(tanpi, a, b))
    @test isnan(enclosemaximum(x -> 1/x, a, b))
end

@testset "discontinuous" begin
    prec = 128
    RR = RealField(prec)
    atol = sqrt(eps())
    rtol = sqrt(eps())

    a = RR(0)
    b = RR(1)

    enclosure = enclosemaximum(a, b, atol = atol, rtol = rtol) do x
        if x < 0.5
            return sin(x)
        elseif x > 0.5
            return cos(x)
        end

        setunion(sin(x), cos(x))
    end
    value = cos(RR(0.5))

    @test contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol
end

@testset "taylor" begin
    prec = 128
    RR = RealField(prec)
    atol = sqrt(eps())
    rtol = sqrt(eps())

    a = RR(0)
    b = RR(1)

    enclosure = enclosemaximum(sin, a, b, atol = atol, rtol = rtol, evaltype = :taylor)
    value = sin(b)
    @test contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol
end
