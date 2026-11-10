using Statistics

@testset "Liquid and vapor arrays" begin

    L, Lmin, ϕ0 = 4.0 + rand(), 0.001, 0.6 + 0.1*rand()
    for closedornot in [true,false]
        X, dXdt, liquid_realratio = randomXp(L,Lmin,closedornot;chargeratio=ϕ0)

        # if this is a closed case, then shift it so that the last
        # slug crosses the line. Ensures more robust testing
        xshift = closedornot ? L - mean(X[end]) : 0.0
        X .= map(u -> (mod(u[1]+xshift,L),mod(u[2]+xshift,L)),X)        

        lslug = mod.(map(u -> u[2],X) - map(u -> u[1],X),L)
        @test all(lslug .> 0)

        lslug2 = OscillatingHeatPipe.XptoLliquidslug(X,L)
        @test lslug == lslug2

        ϕ = sum(lslug)/L
        @test abs(ϕ - ϕ0) < 0.1
        @test ϕ ≈ liquid_realratio

        Xvapor = OscillatingHeatPipe.getXpvapor(X,closedornot)
        lvap = mod.(map(u -> u[2],Xvapor) - map(u -> u[1],Xvapor),L)
        @test all(lvap .> 0)

        lvap2 = OscillatingHeatPipe.XptoLvaporplug(X,L,closedornot)
        
        @test lvap == lvap2

        # Continuity test
        @test sum(lvap) + sum(lslug) ≈ L

        # Test checking positions in liquid slugs.
        np = length(X)

        # Liquid point should return true
        ir = rand(1:np-1)
        @test OscillatingHeatPipe.ifamong(mean(X[ir]),X)

        # Vapor point should return false
        ir = rand(2:np)
        Xi, Xf = X[ir-1][2], X[ir][1]
        @test !OscillatingHeatPipe.ifamong(0.5*(Xi+Xf),X)

        # Test creation of liquid slug arrays
        N = rand(200:300)
        θi = rand()
        Xarray, θarray = OscillatingHeatPipe.constructXarrays(X,N,θi,L)
        # check that all points lie between 0 and L
        @test all(map(u -> all(0 .<= u .<= L),Xarray))
        # both arrays for a slug are the same length
        @test all(map((u,v) -> length(u)==length(v),Xarray,θarray))


        ## test assembly and disassembly 
        
        dXdt .= [(rand(),rand()) for j in 1:np]
        u = OscillatingHeatPipe.X_to_u(X,dXdt)

        ir = rand(1:np)
        Xi, Xf = X[ir]
        dXi, dXf = dXdt[ir]
        @test u[2*ir] == Xf && u[2*ir-1] == Xi
        @test u[2*np + 2*ir] == dXf && u[2*np + 2*ir - 1] == dXi

        M = rand(np)
        Afilm_start = rand(np)
        Afilm_end = rand(np)
        Lfilm_start = rand(np)
        Lfilm_end = rand(np)
        
        # Test that assembly and disassembly are inverses of each other
        u = OscillatingHeatPipe.XMAL_to_u(X,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end)

        X2,dXdt2,M2,Afilm_start2,Afilm_end2,Lfilm_start2,Lfilm_end2 = OscillatingHeatPipe.u_to_XMAL(u)

        @test X2 == X
        @test dXdt2 == dXdt
        @test M2 == M
        @test Afilm_start2 == Afilm_start
        @test Afilm_end2 == Afilm_end
        @test Lfilm_start2 == Lfilm_start
        @test Lfilm_end2 == Lfilm_end

        # test a case in which there is one more X and dXdt element than the others
        # This occurs in an open tube, e.g., vapor plugs at ends and one slug between them
        # so that Xp contains (0,Xp1) and (Xpn,L) to describe the end vapor plugs.
        X = [(rand(),rand()) for j in 1:np+1]
        dXdt = [(rand(),rand()) for j in 1:np+1]
        u = OscillatingHeatPipe.XMAL_to_u(X,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end)
        X2,dXdt2,M2,Afilm_start2,Afilm_end2,Lfilm_start2,Lfilm_end2 = OscillatingHeatPipe.u_to_XMAL(u)

        @test X2 == X
        @test dXdt2 == dXdt
        @test M2 == M
        @test Afilm_start2 == Afilm_start
        @test Afilm_end2 == Afilm_end
        @test Lfilm_start2 == Lfilm_start
        @test Lfilm_end2 == Lfilm_end

        

    end

end

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

@testset "Film area state" begin
    for closed in (false, true), shape in ("square", "circle")
        sys = initialize_ohpsys(sys_plate,p_fluid,power,closedornot=closed,slugnum=3,
                               δfilm_relative=0.04,tubeshape=shape)
        thickness = 0.04 * sys.tube.d / 2
        expected_area = sys.tube.Ac * (1 - (1 - 2thickness / sys.tube.d)^2)
        @test all(sys.vapor.Afilm_start .≈ expected_area)
        @test all(sys.vapor.Afilm_end .≈ expected_area)
        @test !hasproperty(sys.vapor, :δstart)
        @test !hasproperty(sys.vapor, :δend)

        u = newstate(sys)
        mechanical_length = 4length(sys.liquid.Xp) + 5length(sys.vapor.P)
        fields = OscillatingHeatPipe.u_to_XMAL(u[1:mechanical_length])
        @test fields[4] == sys.vapor.Afilm_start
        @test fields[5] == sys.vapor.Afilm_end

        restored = getcurrentsys_nowall!(u,deepcopy(sys))
        @test restored.vapor.Afilm_start == sys.vapor.Afilm_start
        @test restored.vapor.Afilm_end == sys.vapor.Afilm_end
        @test restored.vapor.P ≈ sys.vapor.P
        @test getMtotal(restored) ≈ getMtotal(sys)
        δ_start, δ_end = getδ([restored])
        @test only(δ_start) ≈ thickness
        @test only(δ_end) ≈ thickness
    end
end

@testset "Film areas" begin

    d = sys_tube.tube.d
    Afilm_start = sys_tube.vapor.Afilm_start
    Afilm_end = sys_tube.vapor.Afilm_end
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac

    # Each moving interface below gives a mean slug speed of 0.05.
    Ca = getCa(sys_tube.liquid.μₗ,sys_tube.liquid.σ,0.05)
    δdep = Catoδ(d,Ca,adjust_factor=sys_tube.vapor.ad_fac)

    δarea_start = Afilm_start
    δarea_end = Afilm_end
    δarea_dep = sys_tube.tube.Ac * (1 - (1 - 2δdep / d)^2)
    @test all(Afilm_start .> 0)
    @test all(Afilm_end .> 0)

    # No slugs are moving. Should only be equal to existing films
    Adep = Af_at_interfaces(sys_tube)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Left slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.1,0.0) for i in eachindex(dXdt)]
    Adep = Af_at_interfaces(sys_tube)
    @test all(map(u -> u[1]≈δarea_dep,Adep))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Right slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.0,-0.1) for i in eachindex(dXdt)]
    Adep = Af_at_interfaces(sys_tube)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map(u -> u[2]≈δarea_dep,Adep))

    vol = getVolumevapor(sys_tube)
    @test all(vol .> 0)

    ρv = sys_tube.propconvert.PtoD.(sys_tube.vapor.P)
    @test getMvapor(sys_tube) ≈ ρv.*vol

end


@testset "Mass functions" begin

    numofslugs = length(sys_tube.liquid.Xp)

    sys_tube.vapor.Afilm_start .= sys_tube.tube.peri .* (zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) ) # initial film area
    sys_tube.vapor.Afilm_end .= sys_tube.tube.peri .* (zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) ) # initial film area
    

    d = sys_tube.tube.d
    Afilm_start = sys_tube.vapor.Afilm_start
    Afilm_end = sys_tube.vapor.Afilm_end
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end
    Xp = sys_tube.liquid.Xp
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρₗ
    closedornot = sys_tube.tube.closedornot

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Lliquidslug = XptoLliquidslug(Xp,L)

    ρv = sys_tube.propconvert.PtoD.(sys_tube.vapor.P)

    # mass of vapor
    vol_vapor_analytical = Lvaporplug .* Ac .- Lfilm_start .* Afilm_start .- Lfilm_end .* Afilm_end
    M_vapor_analytical = ρv .* vol_vapor_analytical

    @test OscillatingHeatPipe.getMvapor(sys_tube) ≈ M_vapor_analytical

    # mass of liquid
    vol_liquid_analytical = Lliquidslug .* Ac
    M_liquid_analytical = ρₗ .* vol_liquid_analytical

    @test OscillatingHeatPipe.getMliquid(sys_tube) ≈ M_liquid_analytical

    # mass of films
    vol_film_start_analytical = Lfilm_start .* Afilm_start
    M_film_start_analytical = ρₗ .* vol_film_start_analytical
    vol_film_end_analytical = Lfilm_end .* Afilm_end
    M_film_end_analytical = ρₗ .* vol_film_end_analytical

    Mfilm_start,Mfilm_end = OscillatingHeatPipe.getMfilm(sys_tube)

    @test Mfilm_start ≈ M_film_start_analytical
    @test Mfilm_end ≈ M_film_end_analytical

    @test getMtotal(sys_tube) ≈ sum(M_vapor_analytical) + sum(M_liquid_analytical) + sum(M_film_start_analytical) + sum(M_film_end_analytical)

end

# some threshold values need to be changed for other applications
@testset "Hfilm" begin

    δmin = sys_tube.vapor.δmin;
    δthreshold = 5e-6
    δmax = 1e-4

    kₗ = sys_tube.vapor.k
    Hᵥ = sys_tube.vapor.Hᵥ
    δs = [1e-6;3e-6;1e-5;1.5e-4;3e-4]

    @test OscillatingHeatPipe.Hfilm.(δs,[sys_tube]) ≈ [0.0;(δs[2]-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin);kₗ/δs[3];0.5*kₗ/δmax+1e-6;0.0]
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
                                             # timestep_func=timestep_fourier
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
                    # "angular velocity"         => 0.0,
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
