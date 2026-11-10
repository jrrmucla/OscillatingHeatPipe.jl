using OscillatingHeatPipe
using Random
using Test
using ProgressMeter
using Plots

regression_zero_flux(t, x, base_cache, phys_params, motions) = zeros_surface(base_cache)

function _regression_insulated_polygon(xbound, ybound, Δs)
    body = Polygon(xbound,ybound,Δs)
    
    # stationary rigid-body motion object for the object
    X = MotionTransform([0,0],0)
    joint = Joint(X)
    m = RigidBodyMotion(joint,body)
    x = zero_motion_state(body,m)
    update_body!(body,x,m)

    # no heat flux at the boundaries
    bcdict = Dict("exterior" => regression_zero_flux,"interior" => regression_zero_flux)
    return body, m, bcdict
end

function _get_regression_solid_properties()
    ρₛ = 2730; # material density [kg/m^3]
    cₛ  = 8.93e02; # material specific heat [J/kg K]
    kₛ  = 1.93e02; # material heat conductivity
    αₛ = kₛ/ρₛ/cₛ

    return ρₛ, cₛ, kₛ, αₛ
end

function _place_regression_heaters!(phys_params, power, Δx)
    inches = 2.54e-2
    phys_params["areaheater_power"] = power
    phys_params["areaheater_area"] = 2inches*2inches

    function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
        σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"]
    end

    fregion1_h = Rectangle(0.5inches,1.0inches,1.4*Δx)
    fregion2_h = deepcopy(fregion1_h)
    tr1_h = RigidTransform((0.7inches,-0.0),0.0)
    tr2_h = RigidTransform((-0.7inches,-0.0),0.0)

    update_body!(fregion1_h, tr1_h)
    update_body!(fregion2_h, tr2_h)

    heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!);
    heater2 = AreaForcingModel(fregion2_h,tr2_h,heatermodel!);

    return [heater1, heater2]
end

function _place_regression_coolers!(phys_params, Δx)
    inches = 2.54e-2
    phys_params["areaheater_temp"] = 0.0
    phys_params["areaheater_coeff"] = 2300.0

    function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
        T0 = phys_params["areaheater_temp"]
        h = phys_params["areaheater_coeff"]
        corr = phys_params["flux_correction"]

        σ .= h*(T0 - T) / corr
    end

    fregion1_c = Rectangle(0.5inches,0.0648*0.95/2,1.4*Δx)
    fregion2_c = deepcopy(fregion1_c)
    tr1_c = RigidTransform((2.45inches,-0.0),0.0)
    tr2_c = RigidTransform((-2.45inches,-0.0),0.0)

    update_body!(fregion1_c, tr1_c)
    update_body!(fregion2_c, tr2_c)

    cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!);
    cond2 = AreaForcingModel(fregion2_c,tr2_c,condensermodel!);

    return [cond1, cond2]
end

function _route_regression_ohp(Δx)
    x, y = construct_ohp_curve("ASETS", Δx)
    return BasicBody(x, y)
end

"""
    run_regression_experiment(title, results_dir; seed=1234, inspection=false)

Run the coupled ASETS OHP experiment with repeatable initial slug lengths,
positions, and boiling nucleation sites. Reseed Julia's default RNG immediately
before OHP initialization. Use the same seed when generating reference results.
Save numeric history to `results_dir/title.jld2` and return its path.
Set `inspection=true` to generate GIFs from the live simulation result.
"""
function run_regression_experiment(title, results_dir; seed::Integer=1234, inspection=false)
    Δx=0.0006
    Δt=4E-4
    tf=5.0
    fluid="Butane"

    # 1) Specify Properties
    ρₛ, cₛ, kₛ, αₛ = _get_regression_solid_properties()
    dₛ = 1.5e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)

    ## Fluid Physical Parameters
    Tref = 291.2 # reference temperature, (i.e. freestream temperature)
    p_fluid = SaturationFluidProperty(fluid,Tref) # This function relies on CoolProp.jl package

    # 2) Define Geometry

    ## General Parameters
    Lx = 0.1524; # plate size x [m]
    Ly = 0.0648; # plate size y [m]
    inches=2.54e-2

    ## Simulation Grid
    xlim = (-Lx/2,Lx/2) # plate x limits
    ylim = (-Ly/2,Ly/2) # plate y limits

    g = PhysicalGrid(1.1 .* xlim,1.1 .* ylim,Δx); # build a gird slightly larger than the plate

    ## Conductive Plate

    Δs = 1.4*cellsize(g)

    trim = 0.006
    cond_block = 1.1inches
    
    xbound = [ -Lx/2,-Lx/2, -Lx/2+cond_block,-Lx/2+cond_block, Lx/2-cond_block,  Lx/2-cond_block,
                Lx/2, Lx/2,  Lx/2-cond_block, Lx/2-cond_block,-Lx/2+cond_block, -Lx/2+cond_block];
    ybound = [  Ly/2,-Ly/2, -Ly/2,           -Ly/2+trim,      -Ly/2+trim,       -Ly/2,
            -Ly/2, Ly/2,  Ly/2,            Ly/2-trim ,      Ly/2-trim,        Ly/2 ];
    
    body, m, bcdict = _regression_insulated_polygon(xbound, ybound, Δs)

    ## Evaporators and Condensers
    power = 40 # Watt
    phys_params = Dict( "diffusivity"              => αₛ,
                        "flux_correction"          => ρₛ*cₛ*dₛ,
                        "Fourier"                  => 1.0,
                        "ohp_flux"                 => [NaN], # initial value, the value here is useless
                        "background temperature"   => Tref
                        )

    heaters = _place_regression_heaters!(phys_params, power, Δx)
    coolers = _place_regression_coolers!(phys_params, Δx)

    ## Create Internal OHP Channels
    ohp = _route_regression_ohp(Δx)
    tr_ohp = RigidTransform((0.0,0.0),0.0)

    function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
        σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"]
    end
    
    ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!);

    ## Assemble Forcing Models

    forcing_dict = Dict("heating models" => [heaters...,coolers...,ohp_linesource])

    ## Simulation Timing and Heat-Conduction Problem

    tspan = (0.0, tf); # start time and end time
    dt_record = 0.1   # saving time interval

    timestep_fixed(u,sys) = Δt

    prob = NeumannHeatConductionProblem(
        g,body,scaling=GridScaling,
        phys_params=phys_params,
        bc=bcdict,
        motions=m,
        forcing=forcing_dict,
        timestep_func=timestep_fixed
    );

    ## Construct the Systems

    # Heat-Conduction System
    sys_plate = construct_system(prob)

    # OHP Inner Channel System
    Random.seed!(seed)
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    # Initialize
    u_plate = init_sol(sys_plate)# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate,save_on=false) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

    # Initialize Arrays for Saving
    simulation_result = SimulationResult(integrator_tube,integrator_plate);

    ## Solve
    @showprogress for t in tspan[1]:Δt:tspan[2]
        timemarching!(integrator_tube,integrator_plate,Δt)
        if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
            store!(simulation_result,integrator_tube,integrator_plate)
        end
    end

    # Store simulation data for subsequent regression comparisons.
    mkpath(results_dir)
    save_path = joinpath(results_dir, "$(title).jld2")
    save_results(save_path, simulation_result; kind=:history)
    inspection && write_inspection_files(title, simulation_result, results_dir)

    return save_path
end

function write_inspection_files(title, simulation_result, results_dir)
    t = simulation_result.tube_hist_t

    # Slugs
    anim = @animate for i in eachindex(t)
        plot(OHPSlug(),i,simulation_result)
    end
    gif(anim, joinpath(results_dir, "$(title)_slug.gif"))

    # Temperature
    anim = @animate for i in eachindex(t)
        plot(OHPTemp(),i,simulation_result,clim=(291.2,311.2))
    end
    gif(anim, joinpath(results_dir, "$(title)_temp.gif"))
end

function compare_regression_array(actual, expected; atol, rtol=1e-6)
    @test size(actual) == size(expected)
    size(actual) == size(expected) || return
    @test all(isapprox.(actual, expected; atol=atol, rtol=rtol))
end

"""
    Unpack the state vector and compare each field.
"""
function compare_regression_state(actual, expected)
    # compare separators
    separator = OscillatingHeatPipe.SEPERATION_VAR
    actual_separators = findall(==(separator), actual)
    expected_separators = findall(==(separator), expected)
    @test !isempty(actual_separators)
    @test !isempty(expected_separators)
    @test actual_separators == expected_separators
    isempty(actual_separators) && return
    isempty(expected_separators) && return
    actual_separators == expected_separators || return

    # compare number of entries
    @test length(actual) == length(expected)
    length(actual) == length(expected) || return

    # this exists to support both closed and open tube layouts
    last_mechanical = first(actual_separators) - 1
    @test rem(last_mechanical, 9) in (0, 4)
    rem(last_mechanical, 9) in (0, 4) || return
    a = OscillatingHeatPipe.vectoXMδL(actual[1:last_mechanical])
    b = OscillatingHeatPipe.vectoXMδL(expected[1:last_mechanical])

    # because our values span so many orders of magnitude, 
    # we use different tolerances for each
    names = (:positions, :velocities, :vapor_mass, :film_thickness_start,
             :film_thickness_end, :film_length_start, :film_length_end)
    atols = (1e-9, 1e-8, 1e-12, 1e-11, 1e-11, 1e-9, 1e-9)

    for j in eachindex(names)
        @testset "$(names[j])" begin
            # flatten any tuples we encounter (like interface positions)
            av = j <= 2 ? collect(Iterators.flatten(a[j])) : a[j]
            bv = j <= 2 ? collect(Iterators.flatten(b[j])) : b[j]
            compare_regression_array(av, bv; atol=atols[j])
        end
    end

    @testset "liquid temperatures" begin
        for i in eachindex(actual_separators)
            first_temperature = actual_separators[i] + 1
            last_temperature = i < length(actual_separators) ? actual_separators[i+1] - 1 : length(actual)
            compare_regression_array(actual[first_temperature:last_temperature],
                                     expected[first_temperature:last_temperature]; atol=1e-4, rtol=0)
        end
    end
end

function compare_temperature_values(actual, expected)
    # units: K
    # order of magnitude: ~1E2
    for field in (:tube_hist_θwall, :plate_T_hist)
        @testset "$field" begin
            a = getproperty(actual, field)
            b = getproperty(expected, field)

            @test length(a) == length(actual.tube_hist_t)
            @test length(b) == length(expected.tube_hist_t)
            @test length(a) == length(b)
            length(a) == length(b) || continue
            
            for i in eachindex(a)
                compare_regression_array(a[i], b[i]; atol=1e-4, rtol=0)
            end
        end
    end
end

function compare_state_vectors(actual, expected)
    @testset "state vectors" begin
        a = getproperty(actual, :tube_hist_u)
        b = getproperty(expected, :tube_hist_u)
        @test length(a) == length(actual.tube_hist_t)
        @test length(b) == length(expected.tube_hist_t)
        @test length(a) == length(b)
        length(a) == length(b) || return
        
        for i in eachindex(a)
            compare_regression_state(a[i], b[i])
        end
    end
end

function compare_regression_results(actual, expected)
    # expect number of timesteps and time values to be equal
    @testset "timestamps" begin
        compare_regression_array(actual.tube_hist_t, expected.tube_hist_t;
                                 atol=1e-10, rtol=0)
    end

    compare_temperature_values(actual, expected)
    compare_state_vectors(actual, expected)
end

@testset "regression" begin
    title="regression"

    gold_output_file = normpath(joinpath(@__DIR__, "gold","regression.jld2"))
    results_dir = normpath(joinpath(@__DIR__, "results"))

    # Explicit maintenance mode; normal runs never overwrite the reference.
    # PowerShell: $env:UPDATE_REGRESSION = "true"
    update_reference = get(ENV, "UPDATE_REGRESSION", "false") == "true"
    inspection = get(ENV, "REGRESSION_INSPECTION", "false") == "true"
    test_output_file = run_regression_experiment(title, results_dir; seed=1234, inspection=inspection)
    print(test_output_file)

    test_result = load_results(test_output_file; kind=:history)
    if update_reference
        save_results(gold_output_file, test_result; kind=:history)
        @info "Regenerated numeric regression reference" gold_output_file
    end
    gold_result = load_results(gold_output_file; kind=:history)

    compare_regression_results(test_result, gold_result)
end
