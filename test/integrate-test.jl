@testset "Integration" begin
    res = ArbTools.integrate(AcbField(333), x -> sin(x + exp(x)), 0, 8)
    @test string(res) == "[0.34740017265724780787951215911989312465745625486618018388549271361674821398878532052968510434660 +/- 5.31e-96] + i*0"
end

@testset "Integration tolerance" begin
    C = AcbField(64)
    f(x) = exp(-1000 + x)*sin(10x)
    res = ArbTools.integrate(C, f, 0, 1)
    @test string(res) == "[+/- 1.38e-434] + i*0"
    res = ArbTools.integrate(C, f, 0, 1, abs_tol = 0)
    @test string(res) == "[1.574528586972758e-435 +/- 7.34e-451] + i*0"
    res = ArbTools.integrate(C, f, 0, 1, abs_tol = exp(-1000)/(2.0)^64)
    @test string(res) == "[1.574528586972758e-435 +/- 7.34e-451] + i*0"
    g(x) = exp(1000 + x)*sin(10x)
    res = ArbTools.integrate(C, g, 0, 1, rel_tol = 0)
    @test string(res) == "[+/- 6.94e+433] + i*0"
    res = ArbTools.integrate(C, g, 0, 1, rel_tol = 1e-10)
    @test string(res) == "[6.111029167e+433 +/- 2.66e+423] + i*0"
end

@testset "Integration branch cut" begin
    h(x, analytic = false) = begin
        onemx2 = 1 - x^2
        if analytic && contains_nonpositive(real(onemx2))
            return parent(x)(NaN, NaN)
        else
            return sqrt(onemx2)
        end
    end
    res = ArbTools.integrate(AcbField(64), h, 0, 1)
    @test string(res) == "[0.796113019377273535 +/- 7.17e-19] + i*[+/- 2.32e-31]"
end
