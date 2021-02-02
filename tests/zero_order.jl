@testset "Valuations and indicial equations" begin

    R, x = PolynomialRing(Nemo.QQ, "x")
    g = RationalPS{fmpq}(x, 1 - x)
    f = RationalPS{fmpq}(x, R(1))

    R, v = PolynomialRing(Nemo.QQ, ["x$i" for i in 0:MAX_ORD])
    T, t = PowerSeriesRing(Nemo.QQ, 30, "t"; model=:capped_absolute)

    p = v[2] * (1 - v[1]) - v[1] + v[1]^100

    @test fvaluation(p, f, g) == 100

    p = v[1]^2 - v[2]
    slp = sigma_linear_part(p, f, f^2)
    @test map(x -> get_truncation(x, T), slp) == [2 * t, T(-1)]


    p = v[3] - v[1] * v[2]
    @test Z_linear_part(p, IdentityPS{fmpq}(Nemo.QQ), f + f^2) == -1

    f = AnnihilatorBasedPS{fmpq}(v[2] - v[1]^2 - v[1], f + f^2, [0, 1, 0, 0])
    SD = construct_domain(f)
 
    R, v = PolynomialRing(SD, ["x$i" for i in 0:MAX_ORD]) 
    p = v[2] - v[1] * (gen(SD) + 1)
    @test Z_linear_part(p, ConstPS{fmpq}(Nemo.QQ, 0), f + f^2) == 1
end
