@testset "Difference polynomials arithmetic" begin

   R, v = PolynomialRing(Nemo.QQ, ["x$i" for i in 0:MAX_ORD])

    p = v[1] + 2 * v[2]^2 * v[10]^2

    @test sigma(p, 0) == p
    @test sigma(p) == v[2] + 2 * v[3]^2 * v[11]^2
    @test sigma(p, 2) == v[3] + 2 * v[4]^2 * v[12]^2
    @test leader(p) == v[10]
    @test initial(p) == 2 * v[2]^2
    @test separant(p) == 4 * v[2]^2 * v[10]
    @test rittrank(p) == (9, 2)

    @test reduce(v[3] + 1, v[1] * v[2] - 1) == v[1] + 1

    polys = [
        v[3]^4 + v[1] + 1,
        v[2]^3 * 7 - v[1]^10,
        v[2]^10 * v[4]^3 - 10 * v[3]^3,
        v[4] - 10,
        v[5]^2 * v[4] - v[3]
    ]
    @test select_charset(polys) == [polys[2], polys[4]]
end
