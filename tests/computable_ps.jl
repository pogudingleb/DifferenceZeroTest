@testset "Computable power series" begin

    X = IdentityPS{fmpq}(Nemo.QQ)

    prec = 30
    T, t = PowerSeriesRing(Nemo.QQ, prec, "t"; model=:capped_absolute)

    @test get_truncation(X, T) == t
    @test get_truncation(X / (1 - X), T) == sum([t^i for i in 1:prec])
    @test get_truncation((X * X)(1 + X), T) == 1 + 2 * t + t^2
    
    log = LogPS{fmpq}(Nemo.QQ)

    @test get_truncation(log, T) + O(t^4) == t - t^2 * (1//2) + t^3 * (1//3)

    exp = ExpPS{fmpq}(Nemo.QQ)

    @test get_truncation(log(exp - 1), T) == get_truncation(X, T)

    R, v = PolynomialRing(Nemo.QQ, ["x$i" for i in 0:MAX_ORD])
    X_disguised = AnnihilatorBasedPS{fmpq}(v[2] * (1 - v[1]) - v[1], X / (1 - X), [0, 1, 0, 0])
    @test get_truncation(X_disguised, T) == t

    S_disguised = AnnihilatorBasedPS{fmpq}(v[1] * v[2] + v[1] - v[2], X / (1 - X), [0, 1, 1, 1])
    @test get_truncation(S_disguised, T) == t * inv(1 - t)
end
