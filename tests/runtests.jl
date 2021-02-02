using Test
using TestSetExtensions

using Nemo

include("../src/DifferenceZeroTest.jl")

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

