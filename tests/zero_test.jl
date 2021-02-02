@testset "Basic computation is sigma domains" begin

    FAST_CHECK_ORDER = 0

    X = IdentityPS{fmpq}(Nemo.QQ)
    S = X / (1 - X)

    polyring, (x0, x1, x2) = PolynomialRing(Nemo.QQ, ["x0", "x1", "x2"])
    z = AnnihilatorBasedPS{fmpq}(x1 * (1 - x0) - x0, S, [0, 1, 0, 0])
    R = construct_domain(z)
    z = gen(R)

    @time iszero(sigma(z) * (1 - z) - z)
    @test iszero(sigma(z) * (1 - z) - z)
    @time iszero(sigma(z, 2) * (1 - 2 * z) - z)
    @test iszero(sigma(z, 2) * (1 - 2 * z) - z)

    z2 = AnnihilatorBasedPS{fmpq}(x2 * (1 - 2 * x0) - x0, S, [0, 1, 0, 0])
    R = construct_domain(z2)
    z2 = gen(R)
    @test iszero(sigma(z2) * (1 - z2) - z2)
    @test !iszero(sigma(z2) * (1 - z2) - z2 + z2^2)

    #Rationals = SigmaPSDomain{fmpq}(Nemo.QQ, S, X, x1 * (1 - x0) - x0)
    #RatLog = SigmaPSDomain{SigmaPSDomainElem{fmpq}}(Rationals, LogPS(), R(0))
    #RatLogGamma = SigmaPSDomain{SigmaPSDomainElem{SigmaPSDomainElem{fmpq}}}(RatLog, )
end
