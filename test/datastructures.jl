using Statistics
using Test
using OscillatingHeatPipe

@testset "Liquid and Vapor arrays" begin
    L, Lmin, ϕ0 = 4.0 + rand(), 0.001, 0.6 + 0.1*rand()
    for closedornot in [true,false]
        ζ_L, ζ_R, Q_B, liquid_realratio = randomζ(L,Lmin,closedornot;chargeratio=ϕ0)
        ξ_L, ξ_R = ζ_L, ζ_R; # for uniform area of 1

        # if this is a closed case, then shift it so that the last
        # slug crosses the line. Ensures more robust testing
        xshift = closedornot ? L - (ξ_L[end] + ξ_R[end]) / 2 : 0.0
        ξ_L = mod.(ξ_L .+ xshift, L);
        ξ_R = mod.(ξ_R .+ xshift, L);
        ξ_LR = collect(zip(ξ_L, ξ_R));

        liquid_length = mod.(ξ_R - ξ_L, L)
        @test all(liquid_length .> 0)

        @test liquid_length == ξ_to_liquid_length(ξ_L,ξ_R, L)

        ϕ = sum(liquid_length)/L # actual charge ratio by length
        @test abs(ϕ - ϕ0) < 0.1
        @test ϕ ≈ liquid_realratio
        
        Xvapor = OscillatingHeatPipe.getXpvapor(ξ_LR,closedornot)
        vapor_length = mod.(map(u -> u[2],Xvapor) - map(u -> u[1],Xvapor),L)
        @test all(vapor_length .> 0)
    
        @test vapor_length == ξ_to_vapor_length(ξ_L, ξ_R,L,closedornot)

        # continuity test
        @test sum(vapor_length) + sum(liquid_length) ≈ L

        # Test checking positions in liquid slugs.
        np = length(Q_B)

        # Liquid point should return true
        ir = rand(1:np-1)
        @test OscillatingHeatPipe.ifamong(mean(ξ_LR[ir]),ξ_LR)

        # Vapor point should return false
        ir = rand(2:np)
        Xi, Xf = ξ_LR[ir-1][2], ξ_LR[ir][1]
        @test !OscillatingHeatPipe.ifamong(0.5*(Xi+Xf),ξ_LR)

        # Test creation of liquid slug arrays
        N = rand(200:300)
        θi = rand()
        Xarray, θarray = OscillatingHeatPipe.constructXarrays(ξ_LR,N,θi,L)
        # check that all points lie between 0 and L
        @test all(map(u -> all(0 .<= u .<= L),Xarray))
        # both arrays for a slug are the same length
        @test all(map((u,v) -> length(u)==length(v),Xarray,θarray))
    end
end

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

@testset "Film Area Calculation" begin
    d = sys_tube.tube.Dh
    Af_L = sys_tube.vapor.Af_L
    Af_R = sys_tube.vapor.Af_R
    Q_B = sys_tube.liquid.Q_B
    ξ_L = sys_tube.liquid.ζ_L ./ sys_tube.tube.Ac
    ξ_R = sys_tube.liquid.ζ_R ./ sys_tube.tube.Ac

    # Moving slugs below have a bulk speed of 0.05 m/s.
    Ca = getCa(sys_tube.liquid.μₗ,sys_tube.liquid.σ,0.05)
    δdep = Catoδ(d,Ca,adjust_factor=sys_tube.vapor.ad_fac)

    δarea_start = Af_L
    δarea_end = Af_R
    δarea_dep = sys_tube.tube.Ac * (1 - (1 - 2δdep / d)^2)
    @test all(Af_L .> 0)
    @test all(Af_R .> 0)

    # No slugs are moving. Should only be equal to existing films
    Adep_left, Adep_right = Af_at_interfaces(sys_tube, ξ_L, ξ_R, Q_B ./ sys_tube.tube.Ac)
    @test Adep_left == δarea_end
    @test Adep_right == circshift(δarea_start,-1)

    # Left slug interfaces are advancing. Should be equal to deposited film
    Q_B .= 0.05 * sys_tube.tube.Ac
    Adep_left, Adep_right = Af_at_interfaces(sys_tube, ξ_L, ξ_R, Q_B ./ sys_tube.tube.Ac)
    @test all(x -> x ≈ δarea_dep, Adep_left)
    @test Adep_right == circshift(δarea_start,-1)

    # Right slug interfaces are advancing. Should be equal to deposited film
    Q_B .= -0.05 * sys_tube.tube.Ac
    Adep_left, Adep_right = Af_at_interfaces(sys_tube, ξ_L, ξ_R, Q_B ./ sys_tube.tube.Ac)
    @test Adep_left == δarea_end
    @test all(x -> x ≈ δarea_dep, Adep_right)

    vol = get_vapor_volume(sys_tube)
    @test all(vol .> 0)
end


@testset "Mass Calculation Helpers (Uniform Channel)" begin
    numofslugs = length(sys_tube.liquid.Q_B)

    sys_tube.vapor.Af_L .= sys_tube.tube.Pc .* (zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) ) # initial film area
    sys_tube.vapor.Af_R .= sys_tube.tube.Pc .* (zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) ) # initial film area

    d = sys_tube.tube.Dh
    Af_L = sys_tube.vapor.Af_L
    Af_R = sys_tube.vapor.Af_R
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end
    Xp = OscillatingHeatPipe._ξp(sys_tube)
    Ac = sys_tube.tube.Ac
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρₗ
    closedornot = sys_tube.tube.closedornot

    Lvaporplug = ξ_to_vapor_length(first.(Xp), last.(Xp),L,closedornot)
    Lliquidslug = ξ_to_liquid_length(first.(Xp), last.(Xp), L)

    # mass of vapor (uniform channel)
    # verify that the getMvapor follows this calculation
    dens_vapor = sys_tube.propconvert.PtoD.(sys_tube.vapor.P)
    vol_vapor_analytical = Lvaporplug .* Ac .- Lfilm_start .* Af_L .- Lfilm_end .* Af_R
    mass_vapor_analytical = dens_vapor .* vol_vapor_analytical
    @test OscillatingHeatPipe.getMvapor(sys_tube) ≈ mass_vapor_analytical

    # mass of liquid (uniform channel)
    # verify that the getMliquid follows this calculation
    vol_liquid_analytical = Lliquidslug .* Ac
    M_liquid_analytical = ρₗ .* vol_liquid_analytical
    @test OscillatingHeatPipe.getMliquid(sys_tube) ≈ M_liquid_analytical

    # mass of films (either channel)
    # verify that getMfilm follows this calculation
    vol_film_start_analytical = Lfilm_start .* Af_L
    M_film_start_analytical = ρₗ .* vol_film_start_analytical
    vol_film_end_analytical = Lfilm_end .* Af_R
    M_film_end_analytical = ρₗ .* vol_film_end_analytical

    Mfilm_start, Mfilm_end = OscillatingHeatPipe.getMfilm(sys_tube)
    @test Mfilm_start ≈ M_film_start_analytical
    @test Mfilm_end ≈ M_film_end_analytical

    # verfiy that getMtotal produces the sum of
    # vapor + liquid + film
    @test getMtotal(sys_tube) ≈ sum(mass_vapor_analytical) + sum(M_liquid_analytical) + sum(M_film_start_analytical) + sum(M_film_end_analytical)

end

@testset "Hfilm()" begin
    # 
    δmin = sys_tube.vapor.δmin;
    δthreshold = 5e-6
    δmax = 1e-4

    kₗ = sys_tube.vapor.k
    Hᵥ = sys_tube.vapor.Hᵥ

    # verifies that Hfilm produces values consistent with the manual entires below
    δs = [1e-6;3e-6;1e-5;1.5e-4;3e-4]
    expected_coefficients = [0.0;(δs[2]-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin);kₗ/δs[3];0.5*kₗ/δmax+1e-6;0.0]

    @test OscillatingHeatPipe.Hfilm.(δs,[sys_tube]) ≈ expected_coefficients
end


phys_params["areaheater_coeff"] = 0.0
phys_params["adiabatic_coeff"] = 0.0
phys_params["ohp_flux"] = [0.0]
forcing_dict = Dict("heating models" => [heater1,cond1,cond2,ohp_linesource])

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             timestep_func=timestep_fixed)

sys_plate = construct_system(prob)



@testset "Plate" begin
    ohp = OscillatingHeatPipe._get_ohp_from_forcing_list(sys_plate)
    
    @test ohp == ohp_linesource                                     
end


@testset "Plate IBM heat solver" begin

    tspan = (0.0, 1.0)

    sys_plate = construct_system(prob)
    u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    # #	Single eigenmode decay

    xylim = sys_plate.base_cache.g.xlim
    Δx = sys_plate.base_cache.g.Δx
  
    Lx_test = xylim[1][2] .+ 0.5Δx
    Ly_test = xylim[2][2] .+ 0.5Δx

    # heater case for one step
    T0 = Tref + 0.0 # uniform initial temperature

    # integrator_plate.u .= 0.0 # reset to uniform relative T
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)
    tstep = integrator_plate.p.timestep_func(integrator_plate.u,integrator_plate.p)
    ts = tstep:tstep:tstep

    for i in ts
        step!(integrator_plate, tstep, true)
        push!(T_hist,deepcopy(temperature(integrator_plate)))
    end

    # test mean temperature (energy conservation)

    T_hist_mean = mean(T_hist[end])
    ΔT_mean = T_hist_mean - T0
    ΔT_mean_analytical = power * tstep / (ρₛ * cₛ * 2Lx_test * 2Ly_test * dₛ)

    @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1.5e-2)


# =========================================================================================
    phys_params["areaheater_coeff"] = 4000.0
    phys_params["areaheater_temp"] = -1.0
    sys_plate = construct_system(prob)
    u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    # test T rise/drop for one step with heater and condenser
    T0 = Tref # uniform initial temperature

    integrator_plate.u .= 0.0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)

    ts = tstep:tstep:tstep

    for i in ts
        step!(integrator_plate, tstep, true)
        push!(T_hist,deepcopy(temperature(integrator_plate)))
    end

    maximum_T = maximum(T_hist[end])
    minimum_T = minimum(T_hist[end])
    ΔT_max = maximum_T - T0
    ΔT_min = minimum_T - T0
    total_heater_area = sys_plate.phys_params["areaheater_area"]
    qe = power/total_heater_area
    ΔT_max_analytical = qe * tstep / (ρₛ * cₛ * dₛ)


    hc = sys_plate.phys_params["areaheater_coeff"]
    qc = -hc*(-phys_params["areaheater_temp"])
    ΔT_min_analytical = qc * tstep / (ρₛ * cₛ * dₛ)

    @test isapprox(ΔT_max,ΔT_max_analytical,atol=1e-12,rtol=1e-3)
    @test isapprox(ΔT_min,ΔT_min_analytical,atol=1e-12,rtol=1e-3)

# =========================================================================================
    phys_params["areaheater_coeff"] = 0.0
    phys_params["adiabatic_coeff"] = 400.0
    phys_params["areaheater_temp"] = -1.0

    sys_plate = construct_system(prob)
    u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    # test T rise/drop for one step with heater and condenser
    T0 = Tref # uniform initial temperature

    integrator_plate.u .= 0.0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)

    ts = tstep:tstep:tstep

    for i in ts
        step!(integrator_plate, tstep, true)
        push!(T_hist,deepcopy(temperature(integrator_plate)))
    end

    minimum_T = minimum(T_hist[end])

    ΔT_min = minimum_T - T0

    hc = sys_plate.phys_params["adiabatic_coeff"]
    qc = -hc*(-phys_params["areaheater_temp"])
    ΔT_min_analytical = qc * tstep / (ρₛ * cₛ * dₛ)

    @test isapprox(ΔT_min,ΔT_min_analytical,atol=1e-12,rtol=1e-3)
# =========================================================================================
    q1D_value = 1.0
    L = initialize_ohpsys(sys_plate,p_fluid,power).tube.L

    phys_params["areaheater_power"] = 0.0
    phys_params["areaheater_coeff"] = 0.0
    phys_params["adiabatic_coeff"] = 0.0
    phys_params["areaheater_temp"] = 0.0
    phys_params["ohp_flux"] = [-q1D_value] # W/m^2
    # test energy rise/drop for one step with small heater and line source

    sys_plate = construct_system(prob)
    u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    T0 = Tref + 0.0 # uniform initial temperature


    integrator_plate.u .= 0.0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)

    ts = tstep:tstep:tstep

    for i in ts
        step!(integrator_plate, tstep, true)
        push!(T_hist,deepcopy(temperature(integrator_plate)))
    end

    T_hist_mean = mean(T_hist[end])
    ΔT_mean = T_hist_mean - T0
    ΔT_mean_analytical = (-L*q1D_value) * tstep / (ρₛ * cₛ * 2Lx_test * 2Ly_test * dₛ)

    @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1e-2)
                                                 
end

# reset phys_params
phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*dₛ,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN], # initial value, the value here is useless
                    "areaheater_power"         => power, # total power
                    "areaheater_area"          => areaheater_area, # total area
                    "areaheater_temp"          => 0.0,   # relative temperature compared with "background temperature"
                    "areaheater_coeff"         => hc,
                    "adiabatic_coeff"          => h_adiabatic,
                    "background temperature"   => Tref
                     )

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed)

sys_plate = construct_system(prob)
