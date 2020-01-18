@testset "standard" begin
    prec = 128
    atol = sqrt(eps(ArbReal{prec}))
    rtol = sqrt(eps(ArbReal{prec}))

    a = ArbReal(0, bits = prec)
    b = ArbReal(1, bits = prec)

    enclosure = enclosemaximum(sin, a, b, atol = atol, rtol = rtol)
    value = sin(b)
    @test ArbToolsNemo.contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol

    enclosure = enclosemaximum(sinpi, a, b, atol = atol, rtol = rtol)
    value = ArbReal(1)
    @test ArbToolsNemo.contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol

    enclosure = enclosemaximum(x -> 2x + 1, a, b, atol = atol, rtol = rtol)
    value = 2b + 1
    @test ArbToolsNemo.contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol
end

@testset "unbounded" begin
    prec = 128

    a = ArbReal(0, bits = prec)
    b = ArbReal(1, bits = prec)

    @test isnan(enclosemaximum(tanpi, a, b))
    @test isnan(enclosemaximum(x -> 1/x, a, b))
end

@testset "discontinuous" begin
    prec = 128
    atol = sqrt(eps(ArbReal{prec}))
    rtol = sqrt(eps(ArbReal{prec}))

    a = ArbReal(0, bits = prec)
    b = ArbReal(1, bits = prec)

    enclosure = enclosemaximum(a, b, atol = atol, rtol = rtol) do x
        if x < ArbReal(0.5)
            return sin(x)
        elseif x > ArbReal(0.5)
            return cos(x)
        end

        ArbToolsNemo.setunion(sin(x), cos(x))
    end
    value = sin(ArbReal(0.5))

    @test isnan(enclosemaximum(tanpi, a, b))
    @test isnan(enclosemaximum(x -> 1/x, a, b))
end

@testset "taylor" begin
    prec = 128
    atol = sqrt(eps(ArbReal{prec}))
    rtol = sqrt(eps(ArbReal{prec}))

    a = ArbReal(0, bits = prec)
    b = ArbReal(1, bits = prec)

    sin!(poly, x, n) = begin
        xpoly = ArbToolsNemo._arb_vec_init(2)
        ArbToolsNemo.unsafe_store_ArbRealPtr!(xpoly, x, 1)
        ArbToolsNemo.unsafe_store_ArbRealPtr!(xpoly, ArbReal(1), 2)

        ccall((:_arb_poly_sin_series, :libarb), Cvoid, (Ptr{ArbReal}, Ptr{ArbReal},
                                                        Clong, Clong, Clong),
              poly, xpoly, 2, n, workingprecision(x))

        ArbToolsNemo._arb_vec_clear(xpoly, 2)

        return
    end
    enclosure = enclosemaximum(sin!, a, b, atol = atol, rtol = rtol, evaltype = :taylor)
    value = sin(b)
    @test ArbToolsNemo.contains(enclosure, value)
    @test radius(value) < atol || radius(value)/abs(value) < rtol
end
