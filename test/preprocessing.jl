using Test
using Random
using OscillatingHeatPipe

@testset "Preprocessing: liquid charge volume fraction" begin
    ohp = OscillatingHeatPipe
    L = 1.0
    area = 1e-6
    nslugs = 10
    Random.seed!(1992)

    for closed in (false, true)
        uniform = UniformTube(1e-3, 4e-3, area, L, [0.0, 0.0], closed, 100)
        variable = ohp.VariableTube([0.0, L / 2, L], fill(1e-3, 3),
            fill(4e-3, 3), area .* [1.0, 2.0, 3.0], [0.0, 0.0], closed, 100)

        for (name, tube, volume_to) in (
            ("uniform", uniform, x -> area * x),
            # A(x) = area * (1 + 2x/L), so this integral is independent of the coordinate maps.
            ("variable", variable, x -> area * (x + x^2 / L)))
            total_volume = volume_to(L)
            for charge in (0.25, 0.46, 0.75), sigma in (0.0, 0.01)
                @testset "$name, closed=$closed, charge=$charge, sigma=$sigma" begin
                    left, right, flow, reported = randomζ(tube;
                        numofslugs=nslugs, chargeratio=charge, σ_charge=sigma)
                    ξ_left = ohp.map_ζ_to_ξ(tube, left)
                    ξ_right = ohp.map_ζ_to_ξ(tube, right)
                    volumes = [volume_to(r) - volume_to(l) +
                               (closed && r < l ? total_volume : 0.0)
                               for (l, r) in zip(ξ_left, ξ_right)]
                    actual = sum(volumes) / total_volume
                    @test all(volumes .> 0)
                    @test actual ≈ reported
                    if iszero(sigma)
                        @test actual ≈ charge
                    else
                        # Each sampled volume deviates by at most sqrt(3)*sigma*V/sqrt(N).
                        @test abs(actual - charge) <= sqrt(3 * nslugs) * sigma
                    end
                    @test length(flow) == nslugs
                    @test all(iszero, flow)
                    if !closed
                        @test left[1] == 0
                        @test right[end] ≈ total_volume
                    end
                end
            end
        end
    end
end
