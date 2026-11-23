using Test
using OscillatingHeatPipe
using Random

# A Tale of Two Tubes
Dh = 1E-3
seed = 1234
Q0_B = Dh^2 # QB normally initializes to zero, but that would make continuity trivial thus uninformative here

Random.seed!(seed)
uniform_sys0 = get_adiabatic_system_fixture(Dh, :UNIFORM)
uniform_sys0.liquid.Q_B .= Q0_B

Random.seed!(seed)
variable_sys0 = get_adiabatic_system_fixture(Dh, :VARIABLE)
variable_sys0.liquid.Q_B .= Q0_B

@testset "VariableSystem reproduces UniformSystem when geometry matches" begin
    uniform_u0 = newstate(uniform_sys0)
    variable_u0 = newstate(variable_sys0)

    @testset "newstate()" begin
        @test uniform_u0 ≈ variable_u0
    end

    uniform_du0  = dynamicsmodel(uniform_u0, uniform_sys0)
    variable_du0 = dynamicsmodel(variable_u0, variable_sys0)
    nL = length(uniform_sys0.liquid.Q_B)

    @testset "dynamicsmodel()" begin
        # liquid continuity
        @test uniform_du0[1:nL]     ≈ variable_du0[1:nL]      # Q_left
        @test uniform_du0[nL+1:2nL] ≈ variable_du0[nL+1:2nL]  # Q_right

        # other tests follow as dynamics model is expanded... i.e.
        # liquid momentum
        @test uniform_du0[2nL+1:3nL] ≈ variable_du0[2nL:3nL] broken=true # Q̇_bulk
    end
end