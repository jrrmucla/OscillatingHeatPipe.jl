using Test
using OscillatingHeatPipe

@testset "Circular film geometry" begin
    d = 0.002
    radius = d / 2
    area = π * radius^2
    # Independent annulus geometry, including a film too thick for A ≈ Pδ.
    for fraction in (0.0, 0.04, 0.4, 0.9, 1.0)
        thickness = fraction * radius
        annulus = π * (radius^2 - (radius - thickness)^2)
        @test getδarea(area, d, thickness) ≈ annulus
        @test getδFromδarea(area, d, annulus) ≈ thickness
    end
    tiny = radius * 1e-12
    @test getδFromδarea(area, d, getδarea(area, d, tiny)) ≈ tiny
    @test_throws DomainError getδFromδarea(area, d, 1.01area)

    # Local tube geometry must be used for both initialization and inversion.
    ohp = OscillatingHeatPipe
    diameters = [0.001, 0.002]
    areas = π .* (diameters ./ 2).^2
    tube = ohp.VariableTube([0.0, 1.0], diameters, π .* diameters,
                           areas, [0.0, 0.0], false, 2)
    start, finish = ohp._initial_film_areas(tube, [(0.0, 1.0)], 0.4)
    @test only(start) ≈ areas[1] * (1 - 0.6^2)
    @test only(finish) ≈ areas[2] * (1 - 0.6^2)
end
