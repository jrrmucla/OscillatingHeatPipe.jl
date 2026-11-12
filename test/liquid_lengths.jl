using Test
using OscillatingHeatPipe

@testset "Liquid lengths from separate interface coordinates" begin
    @test ξ_to_liquid_length(1, 4, 10) == 3
    @test ξ_to_liquid_length(9.0, 1.0, 10.0) == 2.0
    @test ξ_to_liquid_length(-1.0, 11.0, 10.0) == 2.0
    @test ξ_to_liquid_length(2.0, 2.0, 10.0) == 0.0
    @test ξ_to_liquid_length(0.0, 10.0, 10.0) == 0.0
    @test ξ_to_liquid_length([1.0, 9.0, 2.0], [4.0, 1.0, 2.0], 10.0) == [3.0, 2.0, 0.0]
    @test isempty(ξ_to_liquid_length(Float64[], Float64[], 10.0))
    @test_throws DimensionMismatch ξ_to_liquid_length([1.0], [2.0, 3.0], 10.0)
end
