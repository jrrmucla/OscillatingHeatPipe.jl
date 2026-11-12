using Test
using Interpolations
using OscillatingHeatPipe
# needs update for new liquid / state vector

@testset "Film interface areas: constant geometry" begin
    ohp = OscillatingHeatPipe
    diameter, perimeter, area = 0.001, 0.004, 1e-6
    positions = [(0.0, 0.1), (0.2, 0.3), (0.4, 0.5), (0.6, 0.7), (0.8, 1.0)]
    placeholder = linear_interpolation([0.0, 1.0], [1.0, 1.0])
    properties = PropConvert("test", placeholder, placeholder, placeholder,
                             placeholder, placeholder)
    wall = Wall(fluid_type="test", boil_type="test", power=0.0,
                L_newbubble=0.01, Xstations=Float64[], boiltime_stations=Float64[],
                Xarray=[0.25, 0.75], θarray=[1.0, 1.0], curvarray=[0.0, 0.0])

    for closed in (false, true)
        nvapor = closed ? length(positions) : length(positions) - 1
        uniform_tube = UniformTube(diameter, perimeter, area, 1.0, [0.0, 0.0], closed, 10)
        variable_tube = ohp.VariableTube([0.0, 0.5, 1.0], fill(diameter, 3),
            fill(perimeter, 3), fill(area, 3), [0.0, 0.0], closed, 10)

        # Cover deposition, existing films, zero velocity, and independent interface signs.
        for (motion, velocities) in (
            ("positive", fill((0.02, 0.02), 5)),
            ("negative", fill((-0.02, -0.02), 5)),
            ("stationary", fill((0.0, 0.0), 5)),
            ("mixed", [(0.02, -0.01), (-0.03, 0.01), (0.0, 0.02),
                       (0.01, 0.0), (-0.02, -0.03)]))
            liquid = Liquid(1.0, 600.0, 2000.0, 1e-7, 1e-4, 0.02,
                copy(positions), copy(velocities), [collect(x) for x in positions],
                [zeros(2) for _ in positions])
            # Distinct film areas expose left/right and vapor-index mistakes.
            vapor = Vapor(k=0.1, P=fill(1e5, nvapor),
                Af_L=collect(1:nvapor) .* 1e-9,
                Af_R=collect(1:nvapor) .* 7e-9,
                Lfilm_start=fill(0.01, nvapor), Lfilm_end=fill(0.01, nvapor))
            uniform_mapping = Mapping(ntuple(_ -> placeholder, 6)...)
            variable_mapping = Mapping(ntuple(_ -> placeholder, 6)...)
            uniform = UniformPhpSystem(uniform_tube, liquid, vapor, wall, properties,
                                       uniform_mapping, Cache(Any[], 0.0))
            variable = ohp.VariableAreaPhpSystem(variable_tube, liquid, vapor, wall,
                properties, variable_mapping, Cache(Any[], 0.0))

            speed = sum(abs((v[1] + v[2]) / 2) for v in velocities) / length(velocities)
            thickness = Catoδ(diameter, getCa(liquid.μₗ, liquid.σ, speed);
                              adjust_factor=vapor.ad_fac)
            deposition = area * (1 - (1 - 2thickness / diameter)^2)
            expected = fill((0.0, 0.0), length(positions))
            active = closed ? eachindex(expected) : (2:length(expected)-1)
            for i in active
                left_vapor = closed ? i : i - 1
                right_vapor = closed ? mod1(i + 1, nvapor) : i
                expected[i] = (
                    velocities[i][1] > 0 ? deposition : vapor.Af_R[left_vapor],
                    velocities[i][2] < 0 ? deposition : vapor.Af_L[right_vapor])
            end
            @testset "closed=$closed, motion=$motion" begin
                for (kind, sys) in (("uniform", uniform), ("variable", variable))
                    @testset "$kind" begin
                        actual = ohp.Af_at_interfaces(sys)
                        @test length(actual) == length(expected)
                        @test all(isapprox.(first.(actual), first.(expected)))
                        @test all(isapprox.(last.(actual), last.(expected)))
                        if !closed
                            @test first(actual) == last(actual) == (0.0, 0.0)
                        end
                    end
                end
            end
        end
    end
end
