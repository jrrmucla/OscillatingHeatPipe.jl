using Test
using Interpolations
using OscillatingHeatPipe

# Exercise mapping assembly without depending on the in-progress fluid structs.
struct GeometryTestSystem{T} <: OscillatingHeatPipe.AbstractPHP
    tube::T
    mapping::Mapping
end
OscillatingHeatPipe.sys_interpolation(sys::GeometryTestSystem) =
    ntuple(_ -> sys.mapping.θ_interp_walltoliquid, 5)

function check_tube_geometry_evaluation(tube::T) where T
    @test (@inferred OscillatingHeatPipe.Dh(tube, 0.5)) isa Float64
    @test (@inferred OscillatingHeatPipe.Pc(tube, 0.5)) isa Float64
    @test (@inferred OscillatingHeatPipe.Ac(tube, 0.5)) isa Float64
    @test (@inferred OscillatingHeatPipe.map_ξ_to_ζ(tube, 0.5)) isa Float64
    @test (@inferred OscillatingHeatPipe.map_ζ_to_ξ(tube, 0.5)) isa Float64
    OscillatingHeatPipe.Dh(tube, 0.5)
    OscillatingHeatPipe.map_ξ_to_ζ(tube, 0.5)
    OscillatingHeatPipe.map_ζ_to_ξ(tube, 0.5)
    @test (@allocated OscillatingHeatPipe.Dh(tube, 0.5)) == 0
    @test (@allocated OscillatingHeatPipe.map_ξ_to_ζ(tube, 0.5)) == 0
    @test (@allocated OscillatingHeatPipe.map_ζ_to_ξ(tube, 0.5)) == 0
end

@testset "Tube-owned geometry" begin
    ohp = OscillatingHeatPipe
    uniform = UniformTube(2.0, 6.0, 3.0, 2.0, [0.0, 0.0], false, 10)
    knots = [0.0, 1.0, 2.0]
    diameters, perimeters, areas = [2.0, 4.0, 2.0], [6.0, 8.0, 6.0], [3.0, 5.0, 3.0]
    variable = ohp.VariableTube(knots, diameters, perimeters, areas,
                                [0.0, 0.0], false, 10)
    @test ohp.Dh(variable, 0.5) ≈ 3.0
    @test ohp.Pc(variable, 0.5) ≈ 7.0
    @test ohp.Ac(variable, 0.5) ≈ 4.0
    @test ohp.map_ξ_to_ζ(variable, 0.5) ≈ 1.75
    @test ohp.construct_geometric_interpolations(variable) ===
        (variable.map_ξ_to_Dh, variable.map_ξ_to_Pc, variable.map_ξ_to_Ac)

    # Mutating caller-owned arrays must not invalidate cached geometry.
    fill!(knots, -1.0)
    fill!(diameters, -1.0)
    fill!(perimeters, -1.0)
    fill!(areas, -1.0)
    @test variable.ξ_geometric_knots == [0.0, 1.0, 2.0]
    @test ohp.Dh(variable, 0.5) ≈ 3.0
    @test ohp.Pc(variable, 0.5) ≈ 7.0
    @test ohp.Ac(variable, 0.5) ≈ 4.0
    @test ohp.map_ξ_to_ζ(variable, 0.5) ≈ 1.75

    placeholder = linear_interpolation([0.0, 2.0], [1.0, 1.0])
    mapping = Mapping(ntuple(_ -> placeholder, 6)...)
    @test fieldcount(Mapping) == 6
    for tube in (uniform, variable)
        check_tube_geometry_evaluation(tube)
        positions = [0.0 0.5; 1.0 2.0]
        for lookup in (ohp.Dh, ohp.Pc, ohp.Ac)
            @test lookup(tube, positions) == lookup.(Ref(tube), positions)
        end
        @test ohp.map_ζ_to_ξ(tube, ohp.map_ξ_to_ζ(tube, positions)) ≈ positions
        forward, inverse = ohp.construct_exact_coordinate_maps(tube)
        @test inverse(forward(positions)) ≈ positions
        sys = GeometryTestSystem(tube, mapping)
        @test ohp.Dh(sys, positions) == ohp.Dh(tube, positions)
        @test ohp.Pc(sys, positions) == ohp.Pc(tube, positions)
        @test ohp.Ac(sys, positions) == ohp.Ac(tube, positions)
        initial = ohp.initialize_interpolations(sys, placeholder)
        refreshed = ohp.refresh_interpolations(sys)
        @test initial isa Mapping
        @test refreshed isa Mapping
        @test refreshed.heightg_interp === placeholder
        @test sys.tube === tube
    end
    @test ohp.Dh(uniform, 0.5) == 2.0
    @test ohp.Pc(uniform, 0.5) == 6.0
    @test ohp.Ac(uniform, 0.5) == 3.0
    @test ohp.map_ξ_to_ζ(uniform, 0.5) == 1.5
end
