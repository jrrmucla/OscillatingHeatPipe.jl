using Test
using Random
using OscillatingHeatPipe

"""
Problem: For non-uniform tubes, the line source in the ILM problem doesn't contain thickness, 
so initialize_ohpsys can't build an accurate tube struct for the system.

Solution: We can give initialize_ohpsys a preexisting tube struct, and it will run verification
to make sure the ILM-measurable parts of the tube are consistent. 

Test: It should succeed if the tube object is consistent with the ILM problem's line source, 
and fail otherwise.
"""

@testset "Supplied tube initialization" begin
    ohp_module = OscillatingHeatPipe
    fluid = SaturationFluidProperty("Butane", 291.2)
    grid = PhysicalGrid((-0.1, 0.1), (-0.1, 0.1), 0.005)
    plate = Rectangle(0.09, 0.09, 0.007)
    zero_flux(t, x, cache, parameters, motions) = zeros_surface(cache)
    boundary_conditions = Dict("exterior" => zero_flux, "interior" => zero_flux)
    line_flux!(flux, temperature, time, cache, parameters) = fill!(flux, 0.0)

    for closed in (false, true)
        x, y, xf, yf, tube = construct_alternating_ohp_curve(
            2, 0.01, 0.1, 0.02, 0.001, -0.03, 0.06, false, false, 0.0,
            0.001, 0.004, 1e-6, 0.0012, 0.005, 1.5e-6; closedornot=closed)
        shape = closed ? BasicBody(x, y) :
            BasicBody([x; xf], [y; yf]; closuretype=parentmodule(BasicBody).OpenBody)
        forcing = LineForcingModel(shape, RigidTransform((0.0, 0.0), 0.0), line_flux!)
        problem = NeumannHeatConductionProblem(grid, plate; scaling=GridScaling,
            phys_params=Dict("diffusivity" => 1e-4, "Fourier" => 1.0), bc=boundary_conditions,
            forcing=Dict("heating models" => [forcing]))
        substrate = construct_system(problem)

        Random.seed!(42)
        uniform = initialize_ohpsys(substrate, fluid, 40.0; closedornot=closed,
                                    slugnum=3, nucleatenum=2)
        @test uniform isa UniformPhpSystem
        @test uniform.tube.N == numpts(shape)
        @test uniform.tube.L ≈ arclength(shape)
        supplied_uniform = UniformTube(0.0012, 0.005, 1.5e-6, arclength(shape),
                                        [1.0, -2.0], closed, numpts(shape))
        supplied = initialize_ohpsys(substrate, fluid, 40.0; tube=supplied_uniform,
                                     slugnum=3, nucleatenum=2)
        @test supplied.tube === supplied_uniform
        @test supplied.liquid.Hₗ ≈ fluid.kₗ / supplied_uniform.d * ohp_module.DEFAULT_NU

        variable = initialize_ohpsys(substrate, fluid, 40.0; tube, slugnum=3, nucleatenum=2)
        @test variable isa ohp_module.VariableAreaPhpSystem
        @test variable.tube === tube
        @test variable.tube.map_ξ_to_Ac !== nothing
        @test variable.cache.mass ≈ getMtotal(variable)
        @test all(isfinite, variable.mapping.H_interp_liquidtowall.(variable.wall.Xarray))
        @test ohp_module.refresh_interpolations(variable) isa ohp_module.Mapping

        wrong_count = UniformTube(0.001, 0.004, 1e-6, arclength(shape),
                                  [0.0, 0.0], closed, numpts(shape) + 1)
        @test_throws ArgumentError initialize_ohpsys(substrate, nothing, 40.0; tube=wrong_count)
        wrong_length = UniformTube(0.001, 0.004, 1e-6, 2arclength(shape),
                                   [0.0, 0.0], closed, numpts(shape))
        @test_throws ArgumentError initialize_ohpsys(substrate, nothing, 40.0; tube=wrong_length)
        shifted_knots = deepcopy(tube)
        shifted_knots.ξ_geometric_knots[2] += 0.0001
        @test_throws ArgumentError initialize_ohpsys(substrate, nothing, 40.0; tube=shifted_knots)
    end

    # An analytic area profile checks exact integration, including the seam.
    tube = ohp_module.VariableTube([0.0, 0.5, 1.0], fill(0.001, 3), fill(0.004, 3),
                                   [1e-6, 2e-6, 1e-6], [0.0, 0.0], true, 10)
    @test ohp_module._tube_region_volume(tube, 0.0, 1.0) ≈ 1.5e-6
    @test ohp_module._tube_region_volume(tube, 0.75, 0.5) ≈ 0.625e-6
end
