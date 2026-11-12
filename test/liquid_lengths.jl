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

@testset "Vapor lengths from separate interface coordinates" begin
    left = [0.0, 0.4, 0.8]
    right = [0.1, 0.6, 1.0]
    @test ξ_to_vapor_length(left, right, 1.0, false) ≈ [0.3, 0.2]
    @test ξ_to_vapor_length(left, right, 1.0, true) ≈ [0.0, 0.3, 0.2]
    @test ξ_to_vapor_length([0.9, 0.4], [0.1, 0.6], 1.0, true) ≈ [0.3, 0.3]
    @test ξ_to_vapor_length([0.9], [0.1], 1.0, true) ≈ [0.8]
    @test ξ_to_vapor_length([0.0, 0.8], [0.2, 1.0], 1.0, false) ≈ [0.6]
    @test_throws DimensionMismatch ξ_to_vapor_length([0.0], [0.1, 0.2], 1.0, true)
    @test_throws ArgumentError ξ_to_vapor_length([0.0], [0.1], 1.0, false)
    @test_throws ArgumentError ξ_to_vapor_length(Float64[], Float64[], 1.0, true)
    @test_throws ArgumentError ξ_to_vapor_length(Float64[], Float64[], 1.0, false)
    for closed in (false, true)
        @test sum(ξ_to_liquid_length(left, right, 1.0)) +
              sum(ξ_to_vapor_length(left, right, 1.0, closed)) ≈ 1.0
    end
end

@testset "Volumes from volume coordinates" begin
    V = 1e-6
    left = V .* [0.0, 0.4, 0.8]
    right = V .* [0.1, 0.6, 1.0]
    @test ζ_to_liquid_volume(0.9V, 0.1V, V) ≈ 0.2V
    @test ζ_to_liquid_volume(left, right, V) ≈ V .* [0.1, 0.2, 0.2]
    @test ζ_to_liquid_volume(0.0, V, V) == 0.0
    @test ζ_to_vapor_volume(left, right, V, false) ≈ V .* [0.3, 0.2]
    @test ζ_to_vapor_volume(left, right, V, true) ≈ V .* [0.0, 0.3, 0.2]
    @test ζ_to_vapor_volume([0.9V], [0.1V], V, true) ≈ [0.8V]
    @test_throws DimensionMismatch ζ_to_liquid_volume([0.0], [0.1V, 0.2V], V)
    @test_throws DimensionMismatch ζ_to_vapor_volume([0.0], [0.1V, 0.2V], V, true)
    @test_throws ArgumentError ζ_to_vapor_volume([0.0], [0.1V], V, false)
    for closed in (false, true)
        @test sum(ζ_to_liquid_volume(left, right, V)) +
              sum(ζ_to_vapor_volume(left, right, V, closed)) ≈ V
    end
end
