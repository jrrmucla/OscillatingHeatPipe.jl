using Test
using OscillatingHeatPipe

@testset "System volumes from ζ coordinates" begin
    ohp = OscillatingHeatPipe
    area = 1e-6
    for closed in (false, true)
        uniform = UniformTube(1e-3, 4e-3, area, 1.0, [0.0, 0.0], closed, 100)
        variable = ohp.VariableTube([0.0, 0.5, 1.0], fill(1e-3, 3),
            fill(4e-3, 3), area .* [1.0, 2.0, 3.0], [0.0, 0.0], closed, 100)
        for (name, tube, integral) in (("uniform", uniform, x -> area*x),
                                      ("variable", variable, x -> area*(x+x^2)))
            @testset "$name, closed=$closed" begin
                left, right = [0.0, 0.4, 0.8], [0.1, 0.6, 1.0]
                liquid = (ζ_L=integral.(left), ζ_R=integral.(right), ρₗ=600.0)
                nv = closed ? 3 : 2
                vapor = (Af_L=fill(1e-8, nv), Af_R=fill(2e-8, nv),
                         Lfilm_start=fill(0.01, nv), Lfilm_end=fill(0.02, nv), P=fill(1e5, nv))
                sys = (tube=tube, liquid=liquid, vapor=vapor, propconvert=(PtoD=p -> 2.0,))
                gaps = [integral(0.4)-integral(0.1), integral(0.8)-integral(0.6)]
                closed && pushfirst!(gaps, 0.0)
                # Avoid a film in the zero-length wraparound gap.
                if closed
                    vapor.Lfilm_start[1] = 0.0
                    vapor.Lfilm_end[1] = 0.0
                end
                films = vapor.Af_L .* vapor.Lfilm_start + vapor.Af_R .* vapor.Lfilm_end
                expected_liquid = integral.(right) - integral.(left)
                @test get_vapor_volume(sys) ≈ gaps - films
                @test getMvapor(sys) ≈ 2 .* (gaps - films)
                @test getMliquid(sys) ≈ 600 .* expected_liquid
                @test sum(get_vapor_volume(sys)) + sum(films) + sum(expected_liquid) ≈ integral(1.0)
                @test ohp._tube_volume(tube) ≈ integral(1.0)
                @test getchargeratio(sys) ≈ getMtotal(sys) / (600 * integral(1.0))
            end
        end
    end
end
