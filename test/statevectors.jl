using Test
using OscillatingHeatPipe

@testset "State Vector Assembly & Disassembly" begin
    for closedornot in (false, true)
        @testset "closed=$closedornot" begin
            tube = UniformTube(1e-3, 4e-3, 1e-6, 1.0, [0.0, 0.0], closedornot, 100)
            ζ_L, ζ_R, Q_B, _ = randomζ(tube;
                numofslugs=4, chargeratio=0.6, σ_charge=0.0)
            nL = length(Q_B)
            nV = closedornot ? nL : nL - 1
            Q_B .= range(-2e-8, 2e-8; length=nL)
            M = collect(1:nV) .* 1e-9
            Af_L = collect(1:nV) .* 2e-8
            Af_R = collect(1:nV) .* 3e-8
            Lfilm_start = collect(1:nV) .* 0.01
            Lfilm_end = collect(1:nV) .* 0.02

            u = OscillatingHeatPipe.ζMAL_to_u(
                ζ_L, ζ_R, Q_B, M, Af_L, Af_R, Lfilm_start, Lfilm_end)

            # expected packing order
            @test length(u) == 3nL + 5nV
            @test u[1:nL] == ζ_L
            @test u[nL+1:2nL] == ζ_R
            @test u[2nL+1:3nL] == Q_B
            @test u[3nL+1:3nL+nV] == M
            @test u[3nL+nV+1:3nL+2nV] == Af_L
            @test u[3nL+2nV+1:3nL+3nV] == Af_R
            @test u[3nL+3nV+1:3nL+4nV] == Lfilm_start
            @test u[3nL+4nV+1:end] == Lfilm_end

            ζ_L2, ζ_R2, Q_B2, M2, Af_L2, Af_R2, Lfilm_start2, Lfilm_end2 =
                OscillatingHeatPipe.u_to_ζMAL(u)

            # inversion
            @test ζ_L2 == ζ_L
            @test ζ_R2 == ζ_R
            @test Q_B2 == Q_B
            @test M2 == M
            @test Af_L2 == Af_L
            @test Af_R2 == Af_R
            @test Lfilm_start2 == Lfilm_start
            @test Lfilm_end2 == Lfilm_end
        end
    end
end
