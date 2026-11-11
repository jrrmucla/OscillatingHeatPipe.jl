using Test
using OscillatingHeatPipe

# Specialize on the callable types so allocation measurements reflect solver
# calls, and warm each operation before measuring it.
function check_coordinate_map_evaluation(forward::F, inverse::G) where {F,G}
    samples = collect(range(0.0, 3.0; length=101))
    @test inverse(forward(samples)) ≈ samples
    @test (@inferred forward(0.5)) isa Float64
    @test (@inferred inverse(0.5)) isa Float64
    forward(0.5)
    inverse(0.5)
    @test (@allocated forward(0.5)) == 0
    @test (@allocated inverse(0.5)) == 0

    dest = similar(samples)
    dest .= forward.(samples)
    @test dest == forward(samples)
    @test (@allocated dest .= forward.(samples)) == 0
end

@testset "Exact coordinate maps" begin
    ohp = OscillatingHeatPipe
    make_tube(knots, areas) = ohp.VariableTube(
        knots, ones(length(knots)), ones(length(knots)),
        areas, [0.0, 0.0], false, length(knots))

    # Constant, increasing, and decreasing area; unequal segment widths.
    for areas in ([2.0, 2.0, 2.0], [1.0, 3.0, 2.0], [3.0, 1.0, 4.0])
        tube = make_tube([0.0, 1.0, 3.0], areas)
        forward, inverse = ohp.construct_exact_coordinate_maps(tube)
        @test forward === tube.map_ξ_to_ζ
        @test inverse === tube.map_ζ_to_ξ
        first_volume = (areas[1] + areas[2]) / 2
        total_volume = first_volume + areas[2] + areas[3]
        @test forward([0.0, 1.0, 3.0]) ≈ [0.0, first_volume, total_volume]
        @test forward(0.5) ≈ areas[1] / 2 + (areas[2] - areas[1]) / 8
        check_coordinate_map_evaluation(forward, inverse)
        @test_throws DomainError forward(-1.0)
        @test_throws DomainError inverse(-1.0)
    end

    @test_throws ArgumentError ohp.construct_exact_coordinate_maps(
        make_tube([0.0, 1.0], [1.0, 0.0]))
end

@testset "Parity between uniform and not-so-variable tubes" begin
    ohp = OscillatingHeatPipe

    Pc = 8.0;
    Ac = 4.0;
    Dh = (4*Ac)/Pc;

    L = 10.0;
    uniform_tube = ohp.UniformTube(Dh, Pc, Ac, L, [0.0, 0.0], false, 2);
    equivalent_variable_tube = ohp.VariableTube(
        [0,L], [Dh,Dh], [Pc,Pc], [Ac,Ac],
        [0.0, 0.0], false, 2);

    ξs = range(0.0, L; length=10)
    ζs = range(0.0, Ac * L; length=10)

    @test ohp.map_ξ_to_ζ(uniform_tube, ξs) ≈ ohp.map_ξ_to_ζ(equivalent_variable_tube, ξs)
    @test ohp.map_ζ_to_ξ(uniform_tube, ζs) ≈ ohp.map_ζ_to_ξ(equivalent_variable_tube, ζs)
end
