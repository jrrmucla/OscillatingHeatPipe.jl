using Test
using Interpolations
using OscillatingHeatPipe

@testset "Shared liquid thermal model" begin
    ohp = OscillatingHeatPipe
    placeholder = linear_interpolation([0.0, 1.0], [10.0, 10.0])
    properties = PropConvert("test", ntuple(_ -> placeholder, 5)...)
    wall = Wall(fluid_type="test", boil_type="test", power=0.0,
        L_newbubble=0.01, Xstations=Float64[], boiltime_stations=Float64[],
        Xarray=[0.0, 1.0], θarray=[10.0, 10.0], curvarray=[0.0, 0.0])
    for closed in (false, true)
        nodes = closed ? [[0.8, 0.9, 0.0, 0.1, 0.2]] : [[0.2, 0.3, 0.4, 0.5, 0.6]]
        temperatures = [[10.0, 10.0, 11.0, 10.0, 10.0]]
        liquid = Liquid(2.0, 3.0, 4.0, 0.01, 1e-4, 0.02,
            [0.0], [0.0], [0.0], nodes, temperatures)
        uniform_tube = UniformTube(1.0, 6.0, 2.0, 1.0, [0.0, 0.0], closed, 10)
        constant_tube = ohp.VariableTube([0.0, 1.0], [1.0, 1.0], [6.0, 6.0],
            [2.0, 2.0], [0.0, 0.0], closed, 10)
        varying_tube = ohp.VariableTube([0.0, 1.0], [1.0, 1.0], [6.0, 8.0],
            [2.0, 5.0], [0.0, 0.0], closed, 10)
        results = Vector{Vector{Float64}}()
        for tube in (uniform_tube, constant_tube, varying_tube)
            ctor = tube isa UniformTube ? UniformPhpSystem : ohp.VariablePhpSystem
            sys = ctor(tube, liquid, Vapor(k=0.1, P=[1e5], Af_L=[0.0], Af_R=[0.0], Lfilm_start=[0.0], Lfilm_end=[0.0]), wall, properties,
                Mapping(ntuple(_ -> placeholder, 6)...), Cache(Any[], 0.0))
            actual = only(liquidmodel(sys))
            # Unit temperature impulse: diffusion into neighbors plus local wall cooling.
            cooling = 2.0 * ohp.Pc(tube, nodes[1][3]) / (3.0 * 4.0 * ohp.Ac(tube, nodes[1][3]))
            @test actual ≈ [0.0, 1.0, -2.0 - cooling, 1.0, 0.0]
            push!(results, actual)
            liquid.θarrays[1] .= 10.0
            @test only(liquidmodel(sys)) ≈ zeros(5)
            liquid.θarrays[1] .= 9.0
            expected_heating = [2.0 * ohp.Pc(tube, x) / (12.0 * ohp.Ac(tube, x)) for x in nodes[1]]
            @test only(liquidmodel(sys)) ≈ expected_heating
            liquid.θarrays[1] .= [10.0, 10.0, 11.0, 10.0, 10.0]
            original_nodes = copy(nodes[1])
            nodes[1][2] += 0.01
            @test_throws ArgumentError liquidmodel(sys)
            nodes[1] .= original_nodes
            nodes[1][2] = nodes[1][1]
            @test_throws ArgumentError liquidmodel(sys)
            nodes[1] .= original_nodes
        end
        @test results[1] ≈ results[2]
    end
end
