export dMdtdynamicsmodel,dynamicsmodel_steadyfilm,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray

function _open_tube_dynamics(du, u, sys, numofliquidslug, Adeposit_left, Adeposit_right,
    Ac, Lliquidslug, Xpvapor, ρₗ, Eratio_plus, Eratio_minus,
    Lvaporplug, Lfilm_start, Lfilm_end, L0threshold_film, L0threshold_pure_vapor,
    Af_L, Af_R, rhs_dXdt, rhs_g, heightg, rhs_press, P)

    #. initialize vapor velocities for left (smaller ξ) and right (larger ξ) side
    numofvaporbubble = numofliquidslug - 1
    v_vapor_left_normal=zeros(numofvaporbubble)
    v_vapor_right_normal=zeros(numofvaporbubble)

    # get average velocities of the liquid slugs
    v_momentum = @view u[2*numofliquidslug+1:2:4*numofliquidslug]
    v_momentum_vapor_end = v_momentum[2:end]
    v_momentum_vapor_start = v_momentum[1:end-1]

    # Continuity of Liquid Slugs to determine interface velocities
    # Ref: Yuxuan 2024 (1)
    v_liquid_left_normal  = v_momentum .+ v_momentum .* Adeposit_left ./ (Ac .- Adeposit_left)
    v_liquid_right_normal = v_momentum .+ v_momentum .* Adeposit_right ./(Ac .- Adeposit_right)
    
    rhs_dLdt = -v_momentum .* (v_liquid_right_normal .- v_liquid_left_normal) ./ Lliquidslug

    # vapor δ
    v_vapor_left_normal = v_liquid_right_normal[1:end-1]
    v_vapor_right_normal = v_liquid_left_normal[2:end]

    Adeposit_vapor_end = Adeposit_left[2:end]
    Adeposit_vapor_start = Adeposit_right[1:end-1]


    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative   = dMdt_latent_end   .- dMdt_latent_end_positive

    dLdt_start,dLdt_end,dAdt_start,dAdt_end,v_vapor_start_final,v_vapor_end_final = film_dynamics(ρₗ, Ac,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
    v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,Adeposit_vapor_start,Adeposit_vapor_end,
    Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Af_L,Af_R)

    v_liquid_left_final = [0.0;vcat(v_vapor_end_final...)]
    v_liquid_right_final = [vcat(v_vapor_start_final...);0.0]

    # up tp here
    for i = 1:numofliquidslug
        du[2*i-1] = v_liquid_left_final[i]
        du[2*i] = v_liquid_right_final[i]

        if (i != 1) && (i != numofliquidslug)
            du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i-1]-P[i]) + rhs_dLdt[i]
        else
            du[2*numofliquidslug + 2*i-1] = 0.0
        end                

        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]
    end


    du[4*numofliquidslug+1:4*numofliquidslug+numofvaporbubble] = dMdt_latent_start+dMdt_latent_end
    du[4*numofliquidslug+numofvaporbubble+1:4*numofliquidslug+2*numofvaporbubble] .= 0.0 # equals to 0 for now
    du[4*numofliquidslug+2*numofvaporbubble+1:4*numofliquidslug+3*numofvaporbubble] .= 0.0 # equals to 0 for now
    # du[4*numofliquidslug+numofvaporbubble+1:4*numofliquidslug+2*numofvaporbubble] = dAdt_start # equals to 0 for now
    # du[4*numofliquidslug+2*numofvaporbubble+1:4*numofliquidslug+3*numofvaporbubble] = dAdt_end # equals to 0 for now
    du[4*numofliquidslug+3*numofvaporbubble+1:4*numofliquidslug+4*numofvaporbubble] = dLdt_start # equals to 0 for now
    du[4*numofliquidslug+4*numofvaporbubble+1:4*numofliquidslug+5*numofvaporbubble] = dLdt_end # equals to 0 for now

    return du
end

function _closed_tube_dynamics(du, u, sys, numofliquidslug, Adeposit_left, Adeposit_right,
    Ac, Lliquidslug, Xpvapor, ρₗ, Eratio_plus, Eratio_minus,
    Lvaporplug, Lfilm_start, Lfilm_end, L0threshold_film, L0threshold_pure_vapor,
    Af_L, Af_R, rhs_dXdt, rhs_g, heightg, rhs_press, P)

    #. initialize vapor velocities for left (smaller ξ) and right (larger ξ) side
    numofvaporbubble = numofliquidslug
    v_vapor_left_normal=zeros(numofvaporbubble)
    v_vapor_right_normal=zeros(numofvaporbubble)

    # get bulk velocities of the liquid slugs
    v_liquid_bulk = @view u[2*numofliquidslug+1:2:4*numofliquidslug]
    v_momentum_vapor_end = v_liquid_bulk
    v_momentum_vapor_start = [v_liquid_bulk[end];v_liquid_bulk[1:end-1]]

    # Continuity of Liquid Slugs to determine interface velocities
    # Ref: Yuxuan 2024 (1)
    v_liquid_left_interface = v_liquid_bulk .+ (v_liquid_bulk .* Adeposit_left ./ (Ac .- Adeposit_left))
    v_liquid_right_interface = v_liquid_bulk .+ (v_liquid_bulk .* Adeposit_right ./(Ac .- Adeposit_right))

    # contraction or extension of the liquid slug from velocity mismatch
    dLdt_liquid = -v_liquid_bulk .* (v_liquid_right_interface .- v_liquid_left_interface) ./ Lliquidslug

    # vapor δ
    v_vapor_left_normal[2:end] = @view v_liquid_right_interface[1:end-1]
    v_vapor_left_normal[1] = v_liquid_right_interface[end]
    v_vapor_right_normal = v_liquid_left_interface

    Adeposit_vapor_end = Adeposit_left
    Adeposit_liquid_end = Adeposit_right
    Adeposit_vapor_start = zeros(size(Adeposit_vapor_end))
    Adeposit_vapor_start[2:end] = @view Adeposit_liquid_end[1:end-1]
    Adeposit_vapor_start[1] = Adeposit_liquid_end[end]

    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative =   dMdt_latent_end   .- dMdt_latent_end_positive

    dLdt_start,dLdt_end,dAdt_start,dAdt_end,v_vapor_start_final,v_vapor_end_final = film_dynamics(ρₗ, Ac,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
    v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,Adeposit_vapor_start,Adeposit_vapor_end,
    Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Af_L,Af_R)

    v_liquid_left_final = v_vapor_end_final
    v_liquid_right_final = [@view v_vapor_start_final[2:end];v_vapor_start_final[1]]

    # this looks like an implementation of the liquid momentum equation (6)
    for i = 1:numofliquidslug
        du[2*i-1] = v_liquid_left_final[i]
        du[2*i] = v_liquid_right_final[i]

        if i != numofliquidslug
            du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + dLdt_liquid[i]
        else
            du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[1]) + dLdt_liquid[i]
        end                

        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

    end

    du[4*numofliquidslug+1:5*numofliquidslug] = dMdt_latent_start+dMdt_latent_end
    du[5*numofliquidslug+1:6*numofliquidslug] = dAdt_start # equals to 0 for now
    du[6*numofliquidslug+1:7*numofliquidslug] = dAdt_end # equals to 0 for now
    du[7*numofliquidslug+1:8*numofliquidslug] = dLdt_start # equals to 0 for now
    du[8*numofliquidslug+1:9*numofliquidslug] = dLdt_end # equals to 0 for now

    return du
end

# this is a equation takes the current u and returns du/dt, p was already updated from the current u.
"""
    This equation takes the current state vector (u) and returns du/dt
"""
function dynamicsmodel(u::Array{Float64,1}, php_sys::UniformPhpSystem)
    
    du = zeros(size(u))

    # extracts essential information stored at sys::UniformPhpSystem
    sys :: UniformPhpSystem = php_sys
    @unpack d,peri,Ac,g,L,closedornot = sys.tube
    @unpack σ,μₗ,ρₗ,Xp,dXdt = sys.liquid
    @unpack P,Eratio_plus,Eratio_minus,Af_L,Af_R,Lfilm_start,Lfilm_end,ad_fac = sys.vapor
    @unpack L_newbubble = sys.wall

    # number of liquid slugs
    numofliquidslug = length(Xp)

    # Yuxuan 2024, equation (3)
    
    Adeposit = Af_at_interfaces(sys)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]


    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Lliquidslug = XptoLliquidslug(Xp,L)

    # get gravity potential
    heightg_interp = sys.mapping.heightg_interp
    Xp1 = [elem[1] for elem in Xp]
    Xp2 = [elem[2] for elem in Xp]
    heightg = map(tuple,heightg_interp(Xp1),heightg_interp(Xp2))

    Xpvapor = getXpvapor(Xp,closedornot)

    # get differential equation factors
    curv_interp = sys.mapping.curv_interp_walltoliquid
    lhs = ρₗ*Ac .* Lliquidslug
    rhs_press = Ac ./ lhs
    V = [mean(elem) for elem in dXdt]
    Re_list = ρₗ .* abs.(V) .* d ./ μₗ
    friction_coefficient = f_churchill.(Re_list .+ 1e-4)
    dXdt_to_stress = -0.125 .* friction_coefficient .* ρₗ .* V .* abs.(V)
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs
    rhs_g = Ac*ρₗ ./ lhs

    # parameters for boolean flags
    L0threshold_film = DEFAULT_L0_FILM_FRAC*L_newbubble
    L0threshold_pure_vapor = DEFAULT_L0_PURE_VAPOR_FRAC*L_newbubble

    if closedornot == false
        return _open_tube_dynamics(du, u, sys, numofliquidslug, Adeposit_left, Adeposit_right,
            Ac, Lliquidslug, Xpvapor, ρₗ, Eratio_plus, Eratio_minus,
            Lvaporplug, Lfilm_start, Lfilm_end, L0threshold_film, L0threshold_pure_vapor,
            Af_L, Af_R, rhs_dXdt, rhs_g, heightg, rhs_press, P)
    end

    if closedornot == true
        return _closed_tube_dynamics(du, u, sys, numofliquidslug, Adeposit_left, Adeposit_right,
            Ac, Lliquidslug, Xpvapor, ρₗ, Eratio_plus, Eratio_minus,
            Lvaporplug, Lfilm_start, Lfilm_end, L0threshold_film, L0threshold_pure_vapor,
            Af_L, Af_R, rhs_dXdt, rhs_g, heightg, rhs_press, P)
    end
end

function film_dynamics(ρₗ, Ac,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
    v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,Adeposit_vapor_start,Adeposit_vapor_end,
    Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Af_L,Af_R)

    # F is the cross-sectional area of the film multiplied by liquid density
    # With area as the state, d(ρ*A*L)/dA = ρ*L.

    F_start = ρₗ .* Af_L
    C_start = ρₗ

    F_end = ρₗ .* Af_R
    C_end = ρₗ

    # the dLdt under normal condition, without considering five different cases
    # first get the mass flow rate for the phase change portion of dLdt, 
    # then convert to the unit of dLdt by F, and then add the interface velocity
    dLdt_start_normal = -(dMdt_latent_start_positive .* Eratio_plus .+ dMdt_latent_start_negative .* Eratio_minus) ./ F_start .- v_vapor_left_normal
    dLdt_end_normal = -(dMdt_latent_end_positive .* Eratio_plus .+ dMdt_latent_end_negative .* Eratio_minus) ./ F_end .+ v_vapor_right_normal

    # get boolean flags for five cases
    he_matrix_start,he_matrix_end = case_flags(Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,dLdt_start_normal,dLdt_end_normal,v_vapor_left_normal,v_vapor_right_normal)

    # println(he_matrix_start)
    # println(he_matrix_end)
    # for case 2, L is very small and the L and δ are frozen to avoid code crash, so the vapor velocity needs to be slightly differ from bulk velocity
    # to account for the mass flow rate caused by the 'frozen' liquid film (very short but still not zero) 's phase change
    v_vapor_left_case2  = v_momentum_vapor_start .+ ((-dMdt_latent_start) ./ (ρₗ .* Ac))
    v_vapor_right_case2 = v_momentum_vapor_end   .- ((-dMdt_latent_end)   ./ (ρₗ .* Ac))

    # for case 5, the liquid interface is essentially eating the other side's liquid film, so need to account for the mass flow rate contributed by that liquid film
    v_vapor_left_case5  = Ac ./ (Ac .- Af_R)   .* v_vapor_left_case2
    v_vapor_right_case5 = Ac ./ (Ac .- Af_L) .* v_vapor_right_case2

    V_vapor_matrix_start = hcat(v_vapor_left_normal,v_vapor_left_case2,v_vapor_left_normal,v_vapor_left_normal,v_vapor_left_case5)'
    V_vapor_matrix_end   = hcat(v_vapor_right_normal,v_vapor_right_case2,v_vapor_right_normal,v_vapor_right_normal,v_vapor_right_case5)'

    v_vapor_start_final = sum(he_matrix_start .* V_vapor_matrix_start,dims=1)
    v_vapor_end_final = sum(he_matrix_end .* V_vapor_matrix_end,dims=1)

    dLdt_matrix_start = hcat(dLdt_start_normal,0 .* dLdt_start_normal,-v_vapor_left_normal,(v_vapor_right_case5 .- v_vapor_left_normal),0 .* dLdt_start_normal)'
    dLdt_matrix_end   = hcat(dLdt_end_normal,0 .* dLdt_end_normal,v_vapor_right_normal,(v_vapor_right_normal .- v_vapor_left_case5),0 .* dLdt_end_normal)'

    dLdt_start = sum(he_matrix_start .* dLdt_matrix_start,dims=1)
    dLdt_end = sum(he_matrix_end .* dLdt_matrix_end,dims=1)

    dAdt_start_normal = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* Adeposit_vapor_start  .* v_vapor_left_normal) ./ (C_start .* Lfilm_start)
    dAdt_end_normal = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* Adeposit_vapor_end .* v_vapor_right_normal) ./ (C_end .* Lfilm_end)

    dAdt_start_case4 = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* Adeposit_vapor_start  .* v_vapor_left_normal .+ ρₗ .* Af_L  .* v_vapor_right_case5) ./ (C_start .* Lfilm_start)
    dAdt_end_case4 = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* Adeposit_vapor_end .* v_vapor_right_normal .- ρₗ .* Af_R  .* v_vapor_left_case5) ./ (C_end .* Lfilm_end)

    dAdt_matrix_start = hcat(dAdt_start_normal,0 .* dAdt_start_normal,dAdt_start_normal,dAdt_start_case4,0 .* dAdt_start_normal)'
    dAdt_matrix_end   = hcat(dAdt_end_normal,0 .* dAdt_end_normal,dAdt_end_normal,dAdt_end_case4,0 .* dLdt_end_normal)'

    dAdt_start = sum(he_matrix_start .* dAdt_matrix_start,dims=1)
    dAdt_end = sum(he_matrix_end .* dAdt_matrix_end,dims=1)

    return dLdt_start,dLdt_end,dAdt_start,dAdt_end,v_vapor_start_final,v_vapor_end_final
end
    
function case_flags(Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,dLdt_start_normal,dLdt_end_normal,v_vapor_left_normal,v_vapor_right_normal)
    # get boolean flags for five cases
    he_evap_start_positive = Bool.(heaviside.(-dLdt_start_normal - v_vapor_left_normal))
    he_evap_end_positive = Bool.(heaviside.(-dLdt_end_normal + v_vapor_right_normal))
    he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
    he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
    he_start_positive = Bool.(heaviside.(dLdt_start_normal))
    he_end_positive = Bool.(heaviside.(dLdt_end_normal))
    he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))


    # zero Lfilm and negative dLdt case
    case2_flag_start = Bool.((1 .- he_meet) .* he_start_short .* (1 .- he_start_positive))
    # two ends meet and both nonzero case
    case3_flag_start = Bool.(he_meet .* he_start_short .* he_start_positive .+ he_meet .* (1 .- he_start_short) .* (1 .- he_end_short)) .* (1 .- he_evap_start_positive) .* (1 .- he_evap_end_positive)
    # two ends meet and other side zero case
    case4_flag_start = Bool.(he_meet .* (1 .- he_start_short) .* he_end_short .* (1 .- he_end_positive))
    # two ends meet and this side zero case
    case5_flag_start = Bool.(he_meet .* he_start_short .* (1 .- he_start_positive))
    # normal case
    case1_flag_start = Bool.((1 .- case3_flag_start) .* (1 .- case2_flag_start) .* (1 .- case4_flag_start) .* (1 .- case5_flag_start))

    # zero Lfilm and negative dLdt case
    case2_flag_end = Bool.((1 .- he_meet) .* he_end_short .* (1 .- he_end_positive))
    # two ends meet and both nonzero case
    case3_flag_end = Bool.(he_meet .* he_end_short .* he_end_positive .+ he_meet .* (1 .- he_end_short) .* (1 .- he_start_short))  .* (1 .- he_evap_start_positive) .* (1 .- he_evap_end_positive)
    # two ends meet and other side zero case
    case4_flag_end = Bool.(he_meet .* (1 .- he_end_short) .* he_start_short .* (1 .- he_start_positive))
    # two ends meet and this side zero case
    case5_flag_end = Bool.(he_meet .* he_end_short .* (1 .- he_end_positive))            
    # normal case
    case1_flag_end = Bool.((1 .- case3_flag_end) .* (1 .- case2_flag_end) .* (1 .- case4_flag_end) .* (1 .- case5_flag_end))

    # he_matrix_start = hcat([case1_flag_start';case2_flag_start';case3_flag_start';case4_flag_start';case5_flag_start'])
    # he_matrix_end   = hcat([case1_flag_end';case2_flag_end';case3_flag_end';case4_flag_end';case5_flag_end'])

    he_matrix_start = hcat(case1_flag_start,case2_flag_start,case3_flag_start,case4_flag_start,case5_flag_start)'
    he_matrix_end   = hcat(case1_flag_end,case2_flag_end,case3_flag_end,case4_flag_end,case5_flag_end)'

    return he_matrix_start,he_matrix_end
end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::UniformPhpSystem)

    dMdt_latent_start=zeros(length(Xpvapor))
    dMdt_latent_end=zeros(length(Xpvapor))
    dMdt_latent_start_positive=zeros(length(Xpvapor))
    dMdt_latent_end_positive=zeros(length(Xpvapor))

    L = sys.tube.L

    peri = sys.tube.Pc
    Ac = sys.tube.Ac
    k = sys.vapor.k

    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    Af_L = sys.vapor.Af_L
    Af_R = sys.vapor.Af_R
    δmin = sys.vapor.δmin

    P = sys.vapor.P
    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    @unpack PtoT,PtoHfg = sys.propconvert
    θ = PtoT.(P)

    Hfg = PtoHfg.(P)

    
    H_interp = sys.mapping.H_interp_liquidtowall
    θ_wall_interp = sys.mapping.θ_interp_walltoliquid

    dx_wall = sys.tube.L/sys.tube.N

    # Lvapor = XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    # Lvapor_pure = max.(Lvapor - Lfilm_start - Lfilm_end,0.0)

    # axialrhs_total = 0

    for i = 1:length(Xpvapor)

        Nstart = Int64(max(2 , div(Lfilm_start[i],dx_wall)))
        
        heatflux_start = quad_trap(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)
        heatflux_start_positive = quad_trap_positive(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))

        heatflux_end = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)
        heatflux_end_positive = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])


        axial_rhs_end = Ac*k*slope_r /Hfg[i]
        axial_rhs_start = Ac*k*(-slope_l) /Hfg[i]

        dMdt_latent_start[i] = heatflux_start*peri/Hfg[i] + axial_rhs_start
        dMdt_latent_end[i] = heatflux_end*peri/Hfg[i] + axial_rhs_end
        dMdt_latent_start_positive[i] = heatflux_start_positive*peri/Hfg[i] + maximum([axial_rhs_start,0.0])
        dMdt_latent_end_positive[i] = heatflux_end_positive*peri/Hfg[i]+ maximum([axial_rhs_end,0.0])

        end

        minimum_film_area = getδarea(sys.tube.Ac, sys.tube.Dh, δmin)
        dMdt_latent_start = dMdt_latent_start.*heaviside.(Af_L .- minimum_film_area)
        dMdt_latent_end = dMdt_latent_end.*heaviside.(Af_R .- minimum_film_area)
        dMdt_latent_start_positive = dMdt_latent_start_positive.*heaviside.(Af_L .- minimum_film_area)
        dMdt_latent_end_positive = dMdt_latent_end_positive.*heaviside.(Af_R .- minimum_film_area)

            return dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive
end

"""
    Models thermal energy within liquid slugs.
    Reference: (9) Yuxuan 2024
"""
function liquidmodel(p::UniformPhpSystem)
    @unpack θarrays,Xarrays,Hₗ,αₗ,Cpₗ,ρₗ = p.liquid
    @unpack Ac,L,peri = p.tube
    @unpack θ_interp_walltoliquid = p.mapping

    du = 0 .* θarrays

    H_rhs = peri / (ρₗ*Cpₗ*Ac)

    for i in eachindex(θarrays)
        
        xs = Xarrays[i];
        dx = mod(xs[end] - xs[1], L) / (length(xs)-1) ## two nodes are at the boundary so need to minus 1 for distance

        fx = map(θ_interp_walltoliquid, xs) - θarrays[i]
        du[i] = αₗ .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs
    end

    return du
end


"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""
function laplacian(u::Vector{Float64})
    unew = zeros(size(u))

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))
    d[1] = -1
    d[end]= -1

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    return (unew)
end


# q'
function sys_to_heatflux(p::UniformPhpSystem)
    @unpack θarray,Xarray = p.wall
    @unpack peri = p.tube
    @unpack θ_interp_liquidtowall,H_interp_liquidtowall = p.mapping

    dθarray = map(θ_interp_liquidtowall, Xarray) .- θarray
    Harray  = map(H_interp_liquidtowall, Xarray)

    return -Harray.*dθarray*peri
end

function integrator_to_heatflux(inte)
    sys_to_heatflux(inte.p)
end

function sys_to_Harray(p::UniformPhpSystem)

    sys = p

    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    Harray  = map(H_interp_liquidtowall, xs)

    Harray
end

function quad_trap(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    L_film_integral = mod((b-a),L)
    h = L_film_integral/N
    mid = mod(a+L_film_integral/2,L)
    H_interp_mid = H_interp(mid)

    int = h * ( H_interp_mid*(θ_interp(a)-θvapor_one) + H_interp_mid*(θ_interp(b)-θvapor_one) ) / 2
    for k=1:N-1
        xk = mod(mod((b-a),L) * k/N + a,L)
        int = int + h*H_interp_mid*(θ_interp(xk)-θvapor_one)
    end
    return int
end

function quad_trap_positive(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    L_film_integral = mod((b-a),L)
    h = L_film_integral/N
    mid = mod(a+L_film_integral/2,L)
    H_interp_mid = H_interp(mid)

    # int = maximum([h * ( H_interp_mid*(θ_interp(a)-θvapor_one) + H_interp_mid*(θ_interp(b)-θvapor_one) ) / 2, 0.0])
    int = maximum([h * ( H_interp_mid*(θ_interp(a)-θvapor_one)/ 2), 0.0]) + maximum([h * ( H_interp_mid*(θ_interp(b)-θvapor_one)/ 2), 0.0])
    for k=1:N-1
        xk = mod(mod((b-a),L) * k/N + a,L)
        int = int + maximum([h*H_interp_mid*(θ_interp(xk)-θvapor_one),0.0])
    end
    return int
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))
