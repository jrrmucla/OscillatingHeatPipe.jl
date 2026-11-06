export run_baseline_experiment, plot_conductance_baseline

include(joinpath(@__DIR__, "helpers.jl"));

using OscillatingHeatPipe # our main package
using Plots               # for plotting
using ProgressMeter       # progress bar in the calculation
using JLD2
import CoolProp

function _get_baseline_solid_properties()
    ρₛ = 2730; # material density [kg/m^3]
    cₛ  = 8.93e02; # material specific heat [J/kg K]
    kₛ  = 1.93e02; # material heat conductivity
    αₛ = kₛ/ρₛ/cₛ

    return ρₛ, cₛ, kₛ, αₛ
end

function _place_baseline_heaters!(phys_params, power, Δx)
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

    heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!);
    heater2 = AreaForcingModel(fregion2_h,tr2_h,heatermodel!);

    return [fregion1_h, fregion2_h], [tr1_h, tr2_h], [heater1, heater2]
end

function _place_baseline_coolers!(phys_params, Δx)
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

    cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!);
    cond2 = AreaForcingModel(fregion2_c,tr2_c,condensermodel!);

    return [fregion1_c, fregion2_c], [tr1_c, tr2_c], [cond1, cond2]
end

function _route_baseline_ohp(Δx)
    x, y = construct_ohp_curve("ASETS", Δx)
    return BasicBody(x, y)
end

function run_baseline_experiment(title; Δx=0.0006, Δt=4E-4, tf=5.0, fluid="Butane")

    # 1) Specify Properties
    ρₛ, cₛ, kₛ, αₛ = _get_baseline_solid_properties()
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
    
    body, m, bcdict = InsulatedStaticPolygon(xbound, ybound, Δs)

    ## Evaporators and Condensers
    power = 40 # Watt
    phys_params = Dict( "diffusivity"              => αₛ,
                        "flux_correction"          => ρₛ*cₛ*dₛ,
                        "Fourier"                  => 1.0,
                        "ohp_flux"                 => [NaN], # initial value, the value here is useless
                        "background temperature"   => Tref
                        )

    fregions_h, trs_h, heaters = _place_baseline_heaters!(phys_params, power, Δx)
    fregions_c, trs_c, coolers = _place_baseline_coolers!(phys_params, Δx)

    ## Create Internal OHP Channels
    
    ohp = _route_baseline_ohp(Δx)
    tr_ohp = RigidTransform((0.0,0.0),0.0)

    function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
        σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"]
    end
    
    ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!);

    ## Inspect Layout

    println("Plotting layout for inspection...")

    plot(body,fillalpha=0)
    for (region, transform) in zip(fregions_h, trs_h)
        update_body!(region, transform)
        plot!(region,fillcolor=:red,fillalpha=0.5)
    end
    for (region, transform) in zip(fregions_c, trs_c)
        update_body!(region, transform)
        plot!(region,fillcolor=:blue,fillalpha=0.5)
    end

    plot!(ohp,closedornot=true,fillalpha=0,color=:black)
    results_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_dir)
    savefig(joinpath(results_dir, "$(title)_layout.png"))

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
    println("Constructing simulation system...")

    # Heat-Conduction System
    sys_plate = construct_system(prob)

    # OHP Inner Channel System
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    # Initialize
    u_plate = init_sol(sys_plate)# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate,save_on=false) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

    # Initialize Arrays for Saving
    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    ## Solve
    @showprogress for t in tspan[1]:Δt:tspan[2]
        timemarching!(integrator_tube,integrator_plate,Δt)
        if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
            store!(SimuResult,integrator_tube,integrator_plate)
        end
    end

    # Store Data
    save_path = joinpath(results_dir, "$(title).jld2")
    save(save_path,"SimulationResult",SimuResult)

    return save_path
end

function plot_conductance_baseline(results::AbstractVector, titles::AbstractVector)
    length(results) == length(titles) || throw(ArgumentError("Provide one title per result."))
    isempty(results) && throw(ArgumentError("Provide at least one result."))

    # Virtual plate sensors at the ASETS-II locations from the post-processing example.
    x2Dsensors = [-2.75,-1.4,-0.8,0.0,0.0,0.8,1.4,2.75] .* INCHES
    y2Dsensors = [0.0,   0.0, 0.0,0.0,0.4,0.0,0.0,0.0] .* INCHES
    plate_sensors = (x2Dsensors,y2Dsensors)
    ihot = 4
    icold = 8
    plt = plot(title="Thermal conductance", xlabel="Time [s]",
               ylabel="G [W/K]", legend=:best)
    for (result, title) in zip(results, titles)
        t_hist,g_hist = getTcurve(plate_sensors,result)
        power = result.integrator_tube.p.wall.power
        keep = t_hist .>= 2.0
        conductance = power ./ (g_hist[keep,ihot] .- g_hist[keep,icold])
        plot!(plt,t_hist[keep],conductance,label=string(title),linewidth=2)
    end
    return plt
end

plot_conductance_baseline(result::SimulationResult, title="Baseline") =
    plot_conductance_baseline([result], [title])
