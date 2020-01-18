@testset "isolateroots-ball" begin
    f = cos
    a = ArbReal(0)
    b = ArbReal(10)
    roots = [ArbReal(π)/2, 3ArbReal(π)/2, 5ArbReal(π)/2]
    found, flags = isolateroots(f, a, b, evaltype = :ball)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test !any(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))

    f = x -> (x - 1)*(x - 2)*(x - 3)
    a = ArbReal(0)
    b = ArbReal(5)
    roots = ArbReal.([1, 2, 3])

    found, flags = isolateroots(f, a, b, evaltype = :ball)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test !any(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))

    f = x -> exp(x) - 1
    a = ArbReal(0)
    b = ArbReal(5)
    roots = ArbReal.([0])

    found, flags = isolateroots(f, a, b, evaltype = :ball)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test !any(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))
end

@testset "isolateroots-taylor" begin
    f! = (poly, x, n) -> begin
        ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, cos(x), 1)
        if n > 1
            ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, -sin(x), 2)
        end

        nothing
    end
    a = ArbReal(0)
    b = ArbReal(10)
    roots = [ArbReal(π)/2, 3ArbReal(π)/2, 5ArbReal(π)/2]
    found, flags = isolateroots(f!, a, b, evaltype = :taylor)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test all(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))

    f! = (poly, x, n) -> begin
        ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, (x - 1)*(x - 2)*(x - 3), 1)
        if n > 1
            ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, 3x^2 - 12x + 11, 2)
        end

        nothing
    end
    a = ArbReal(0)
    b = ArbReal(5)
    roots = ArbReal.([1, 2, 3])

    found, flags = isolateroots(f!, a, b, evaltype = :taylor)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test all(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))

    f! = (poly, x, n) -> begin
        ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, exp(x) - 1, 1)
        if n > 1
            ArbToolsNemo.unsafe_store_ArbRealPtr!(poly, exp(x), 2)
        end

        nothing
    end
    a = ArbReal(0)
    b = ArbReal(5)
    roots = ArbReal.([0])

    found, flags = isolateroots(f!, a, b, evaltype = :taylor)

    @test length(found) == length(roots) && length(flags) == length(roots)
    @test all(flags)
    @test all(ArbToolsNemo.contains(setinterval(found[i]...), roots[i]) for i in 1:length(roots))
end
