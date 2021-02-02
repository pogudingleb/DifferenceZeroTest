@testset "Basic computation is sigma domains" begin

    R, x = PolynomialRing(Nemo.QQ, "x")
    g = RationalPS{fmpq}(x, 1 - x)

    S, (x0, x1) = PolynomialRing(Nemo.QQ, ["x0", "x1"])
    X = AnnihilatorBasedPS{fmpq}(x1 * (1 - x0) - x0, g, [0, 1, 0, 0])
    R2 = construct_domain(X)

    prec = 30
    T, t = PowerSeriesRing(Nemo.QQ, prec, "t"; model=:capped_absolute)

    f = gen(R2)
    @test get_truncation(f, T) == t
    @test get_truncation(sigma(f), T) == sum([t^i for i in 1:prec])
    @test get_truncation(sigma(f) * (1 - f) - f, T) == 0

    R2 = construct_domain(g, LogPS{fmpq}(Nemo.QQ))
    f = gen(R2)

    @test get_truncation(f, T) == sum([t^i * ((-1)^(i - 1) // i) for i in 1:prec])

end
