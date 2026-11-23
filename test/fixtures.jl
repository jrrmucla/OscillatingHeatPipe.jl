using Test
using OscillatingHeatPipe

@testset "Plate Fixtures" begin
    heated = get_powered_plate_fixture()
    insulated = get_adiabatic_plate_fixture()
    @test heated.phys_params["areaheater_power"] == 30.0
    @test insulated.phys_params["areaheater_power"] == 0.0
    overrides = Dict("ohp_flux" => [0.0])
    custom = _construct_plate_fixture(; phys_params=overrides)
    custom.phys_params["ohp_flux"][1] = 1.0
    @test overrides["ohp_flux"] == [0.0]
end

@testset "System Fixtures" begin
    uniform = get_adiabatic_system_fixture(1E-3, :UNIFORM; slugnum=3)
    variable = get_adiabatic_system_fixture(1E-3, :VARIABLE; slugnum=3)
    @test uniform isa OscillatingHeatPipe.UniformPhpSystem
    @test variable isa OscillatingHeatPipe.VariablePhpSystem
    @test variable.tube.Ac_geometric_knots == fill(uniform.tube.Ac, 2)
    @test variable.tube.Dh_geometric_knots == fill(uniform.tube.Dh, 2)
    @test variable.tube.Pc_geometric_knots == fill(uniform.tube.Pc, 2)
    @test last(variable.tube.ξ_geometric_knots) == uniform.tube.L
    @test variable.tube.N == uniform.tube.N
end

sys_plate = get_powered_plate_fixture()
phys_params = deepcopy(sys_plate.phys_params)
power = phys_params["areaheater_power"]
areaheater_area = phys_params["areaheater_area"]
