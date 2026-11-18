#=
boiling.jl

The OHP model needs a way to handle the formation of new vapor bubbles, this 
callback is responsible for the nitty-gritty details of where, why, and how
those new vapor bubbles form.

We allow boiling only at particular stations that are created during 
initialization, defined by
- Position in ξ: Xstations
- Last time that the station was in vapor: boiltime_stations   
    (including as a result of boiling at said station)

When a new bubble is created...
- length is determined by L_newbubble, a constant in sys.wall
- temperature is taken either as wall temp, or liquid temp -> determines pressure
- film characteristics follow some preset rules, TODO definitely should talk to Yuxuan about
  the basic justification behind them
=#

export boiling_affect!,nucleateboiling,boiling_condition

"""   
Determines whether the integrator should check for boiling events
"""
function boiling_condition(u,t,integrator)
    _reject_nonuniform(integrator.p)
    t_interval = DEFAULT_BOIL_SCAN_INTERVAL
    ϵ = 1e-5
    return (abs(mod(t,t_interval)-t_interval) < ϵ) || mod(t,t_interval) < ϵ
end

"""
Entry point for the boiling callback to affect a PHPSystem integrator
"""
function boiling_affect!(integrator)

    sys = deepcopy(getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p)))

    @unpack wall, propconvert, tube, vapor, liquid = sys
    @unpack fluid_type, PtoT, TtoP, PtoD = propconvert
    @unpack boil_type, boil_interval, Rn, Xstations, boiltime_stations, L_newbubble = wall
    @unpack Dh = tube
    @unpack P = vapor

    Tavg = (PtoT.(maximum(P)) + PtoT.(minimum(P)))/2

    Δθthreshold = RntoΔT(Rn,Tavg,fluid_type,Dh,TtoP)

    Δθ_array = getsuperheat.(Xstations,Ref(sys))
    superheat_flag = (Δθ_array .> Δθthreshold) .* ((integrator.t .- sys.wall.boiltime_stations) .> boil_interval)

    boil_event_count = 0;
    boiltime_update_flags = Bool.(zero(boiltime_stations))

    # requirements for boiling to occur
    superheat_at_station(i::Int) = superheat_flag[i]
    liquid_at_station(i::Int) = ifamong(Xstations[i], _ξp(sys))
    space_for_new_bubble(i::Int) = is_there_space_for_new_bubble(sys,i)

    # in the loop below, _ξp(sys) cannot be cached, because values will change within the loop
    for i in eachindex(Xstations)
        if liquid_at_station(i) && space_for_new_bubble(i) && superheat_at_station(i)
            push!(integrator.p.cache.boil_hist,[i,integrator.t]);
            boil_event_count += 1;

            if boil_type == "liquid T"
                new_vapor_pressure = sys.mapping.P_interp_liquidtowall(Xstations[i])
            elseif boil_type == "wall T"
                new_vapor_pressure = TtoP(sys.mapping.θ_interp_walltoliquid(Xstations[i]))
            end

            # TODO: ask Yuxuan about this comment: 'P needs to be given from energy equation'
            # println("pre-boiling energy",get_vapor_energy(p)+get_liquid_energy(p))
            sys = nucleateboiling(sys,(Xstations[i]-L_newbubble/2,Xstations[i]+L_newbubble/2),new_vapor_pressure)
            # println("post-boiling energy:",get_vapor_energy(p)+get_liquid_energy(p))

            boiltime_update_flags[i] = true

        elseif !liquid_at_station(i)
            boiltime_update_flags[i] = true
        end
    end

    boiltime_stations_new = boiltime_stations + boiltime_update_flags .* (integrator.t .- boiltime_stations)
    integrator.p.wall.boiltime_stations = boiltime_stations_new

    postboil_vapor_mass = getMvapor(sys)

    unew=[ζMAL_to_u(sys.liquid.ζ_L, sys.liquid.ζ_R, sys.liquid.Q_B, 
        postboil_vapor_mass, sys.vapor.Af_L, sys.vapor.Af_R, sys.vapor.Lfilm_start, sys.vapor.Lfilm_end);
        liquidθtovec(sys.liquid.θarrays)];

    resize!(integrator.u,length(unew))
    integrator.u = deepcopy(unew)
end

function nucleateboiling(sys::UniformPhpSystem, new_vapor_ξ::Tuple{<:Real,<:Real}, new_vapor_P::Real)
    # TODO: understand (and hopefully improve on) the approach here, in interest of time I am leaving as much as possible as-is
    ρₗ = deepcopy(sys.liquid.ρₗ)
    Ac = sys.tube.Ac
    ξp = _ξp(sys)
    Q_B = copy(sys.liquid.Q_B)
    Af_L = deepcopy(sys.vapor.Af_L)
    Af_R = deepcopy(sys.vapor.Af_R)

    P = deepcopy(sys.vapor.P)
    L = sys.tube.L

    ξarrays = sys.liquid.ξarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot
    Lfilm_start = deepcopy(sys.vapor.Lfilm_start)
    Lfilm_end = deepcopy(sys.vapor.Lfilm_end)


    index = getinsertindex(ξp,(new_vapor_ξ[2]+new_vapor_ξ[1])/2,sys.tube.L,sys.tube.closedornot)

    Linsert = mod(new_vapor_ξ[end] - new_vapor_ξ[1],L)

    """let's do constant film thickness for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    δdeposit = 2e-5
    Af_L_new = insert!(Af_L,index+1,getδarea(sys.tube.Ac, sys.tube.Dh, δdeposit))
    Af_R_new = insert!(Af_R,index+1,getδarea(sys.tube.Ac, sys.tube.Dh, δdeposit))


    Nvapor = length(P)
    loop_minus_index_new = [Nvapor+1;1:Nvapor]

    Lfilm_start_new = insert!(Lfilm_start,index+1,Linsert/8)
    Lfilm_end_new = insert!(Lfilm_end,index+1,Linsert/8)

    Lliquid_adjust = 0
    Xpnew = getnewXp(ξp,index,new_vapor_ξ,Lliquid_adjust,L,closedornot)

    # const P in adjacent vapors
    Pnew = insert!(P,index+1,new_vapor_P)

    Q_B_new = copy(Q_B) # momentum conservation
    insert!(Q_B_new,index+1,Q_B_new[index])

    sysnew = deepcopy(sys)

    _set_ξp!(sysnew, Xpnew)
    sysnew.liquid.Q_B = Q_B_new
    sysnew.vapor.P = Pnew
    sysnew.vapor.Af_L = Af_L_new
    sysnew.vapor.Af_R = Af_R_new
    sysnew.vapor.Lfilm_start = Lfilm_start_new
    sysnew.vapor.Lfilm_end = Lfilm_end_new

    # multiply a factor to adjust the length of vapor/liquid to conserve mass, accounting for the inserted vapor volume.
    PtoD = sys.propconvert.PtoD
    ρinsert = PtoD(new_vapor_P)
    Mvapor_new = sum(getMvapor(sysnew))
    Mfilm_new = sum(sum.(getMfilm(sysnew)))
    Mliquid_new = sum(getMliquid(sysnew))
    Mold = sysnew.cache.mass
    L_adjust = (Mvapor_new + Mfilm_new + Mliquid_new - Mold)./ (ρₗ*Ac) .* (1 +  ρinsert/ρₗ)

    Lvaporplug = ξ_to_vapor_length(first.(Xpnew), last.(Xpnew),sysnew.tube.L,sysnew.tube.closedornot)
    Lpurevapor = Lvaporplug .- Lfilm_start_new .- Lfilm_end_new
    Lliquidslug = ξ_to_liquid_length(first.(Xpnew), last.(Xpnew), sys.tube.L)


    maxindex = 0;
    maxvalue_type = 0;
    # L_adjust < 0 # means we need to shrink the vapor part, I sill first try to shrink the pure vapor part, then the film part.
    if L_adjust < 0

        # This part is to limit the shrinking length, so that it will not shrink too much and look unphysical for one boiling event.
        L_adjust_max = DEFAULT_L_ADJUST_MAX_FACTOR*sysnew.wall.L_newbubble
        L_adjust = (-L_adjust > L_adjust_max) ? -L_adjust_max : L_adjust

        # This part try to find the maximum length of the pure vapor, and try to shrink it first.
        Lpurevapor_maxvalueindex = findmax(Lpurevapor)
        Lpurevapor_maxvalue = Lpurevapor_maxvalueindex[1]
        Lpurevapor_maxindex = Lpurevapor_maxvalueindex[2]

        # If the longest pure vapor part is long enough we will shrink it first.
        if Lpurevapor_maxvalue > -L_adjust
            maxvalue = Lpurevapor_maxvalue
            maxindex = Lpurevapor_maxindex
            maxvalue_type = 3 # 1: Lfilm_start, 2: Lfilm_end , 3: Lpurevapor

            Xpnew[maxindex] = mod.((Xpnew[maxindex][1]+L_adjust,Xpnew[maxindex][2]),L)
        # If the longest pure vapor part is not long enough we will shrink the longest liquid film part.
        else
            # find the index of the longest total liquid film part, then find the side of the longest liquid film part, then shrink it.
            maxindex = findmax(Lfilm_start_new .+ Lfilm_end_new)[2]

            Lfilm_maxvalueindex = findmax([Lfilm_start_new[maxindex],Lfilm_end_new[maxindex]])
            maxvalue = Lfilm_maxvalueindex[1]
            maxvalue_type = Lfilm_maxvalueindex[2] # find the section with largest length for future shrinking, 1: Lfilm_start, 2: Lfilm_end, 3: Lpurevapor

            if maxvalue_type == 1
                Afilm = Af_L_new[maxindex]
                factor_adjust = Ac/(Ac-Afilm)

                if maxvalue > -L_adjust*factor_adjust

                sysnew.vapor.Lfilm_start[maxindex] = sysnew.vapor.Lfilm_start[maxindex] + L_adjust*factor_adjust
                Xpnew[loop_minus_index_new[maxindex]] = mod.((Xpnew[loop_minus_index_new[maxindex]][1],Xpnew[loop_minus_index_new[maxindex]][2]-L_adjust*factor_adjust),L)

                else
                    println("L_adjust:",L_adjust)
                    println("maxvalue:",maxvalue)
                    println("maxvalue_type:",maxvalue_type)
                    println("factor_adjust:",factor_adjust)
                    println("case1: mass conservation enforcement failed in this boiling event! if this happens ocassionally, it is fine, but if it happens frequently, your mass conservation may be off so check it!")
                    throw("boiling function cannot shrink or expand a liquid slug")
                end

            elseif maxvalue_type == 2
                Afilm = Af_R_new[maxindex]
                factor_adjust = Ac/(Ac-Afilm)

                if maxvalue > -L_adjust*factor_adjust

                sysnew.vapor.Lfilm_end[maxindex] = sysnew.vapor.Lfilm_end[maxindex] + L_adjust*factor_adjust
                Xpnew[maxindex] = mod.((Xpnew[maxindex][1]+L_adjust*factor_adjust,Xpnew[maxindex][2]),L)

                else
                    println("L_adjust:",L_adjust)
                    println("maxvalue:",maxvalue)
                    println("maxvalue_type:",maxvalue_type)
                    println("factor_adjust:",factor_adjust)
                    println("case2: mass conservation enforcement failed in this boiling event! if this happens ocassionally, it is fine, but if it happens frequently, your mass conservation may be off so check it!")
                    throw("boiling function cannot shrink or expand a liquid slug")
                end
            else
                throw("error in liquid merging, maxvalue_type is not 1,2,3")
            end
        end


        # L_adjust > 0 # means we need to increase the vapor part, I will always increase the pure vapor part.
    else
        # L_newbubble = sysnew.wall.L_newbubble
        # L_adjust = (L_adjust > L_newbubble) ? L_newbubble : L_adjust

        Lliquid_maxvalueindex = findmax(Lliquidslug)
        maxvalue = Lliquid_maxvalueindex[1]
        maxindex = Lliquid_maxvalueindex[2]
        maxvalue_type = 1 # when vapor part needs to increase length, we only increase the Lpurevapor, 1: Lpurevapor, 2: Lfilm_start, 3: Lfilm_end

        if maxvalue > L_adjust
            Xpnew[maxindex] = mod.((Xpnew[maxindex][1]+L_adjust,Xpnew[maxindex][2]),L)
        else
            # maxindex = 0
            println("L_adjust:",L_adjust)
            println("maxvalue:",maxvalue)
            println("maxvalue_type:",maxvalue_type)
            println("factor_adjust:",factor_adjust)
            println("case3: mass conservation enforcement failed in this boiling event! if this happens ocassionally, it is fine, but if it happens frequently, your mass conservation may be off so check it!")
            throw("boiling function cannot shrink or expand a liquid slug")
        end
    end

    # up to now we got the correct Xpnew, next step is to get Xarraysnew, the splitted Xarrays.

    _set_ξp!(sysnew, Xpnew)
    Xarraysnew,θarraysnew = getnewXθarrays(index,Xpnew,ξarrays,θarrays,L,maxindex,maxvalue_type)
    # θarraysnew = getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)

    sysnew.liquid.ξarrays = Xarraysnew
    sysnew.liquid.θarrays = θarraysnew

    sysnew.mapping = refresh_interpolations(sysnew)

    return sysnew
end

function getnewXθarrays(index,Xpnew,Xarrays_old,θarrays_old,L,maxindex,maxvalue_type)
    Nvaporplug = length(Xpnew)
    maxindex_minus = (maxindex != 1) ? maxindex-1 : Nvaporplug

    Nold= length(Xarrays_old[index])

    Xarraysnew = deepcopy(Xarrays_old)
    θarraysnew = deepcopy(θarrays_old)

    Lliquidslug_new_left =  mod((Xpnew[index][end] - Xpnew[index][1]),L)
    Lliquidslug_new_right = mod((Xpnew[index+1][end] - Xpnew[index+1][1]),L)

    arrayindex = ceil.(Int, Nold * Lliquidslug_new_left / (Lliquidslug_new_left + Lliquidslug_new_right))
    arrayindex = (arrayindex != 1) ? arrayindex : 2
    arrayindex = (arrayindex != Nold) ? arrayindex : Nold-1

    Xarraysnewleft = constructoneXarray(Xpnew[index],arrayindex,L)
    Xarraysnewright = constructoneXarray(Xpnew[index+1],Nold-arrayindex+1,L)

    θarraysnewleft = θarrays_old[index][1:arrayindex]
    θarraysnewright= θarrays_old[index][arrayindex+1:end]
    insert!(θarraysnewright, 1, θarrays_old[index][arrayindex])

    splice!(Xarraysnew, index)
    insert!(Xarraysnew, index,Xarraysnewleft)
    insert!(Xarraysnew, index+1,Xarraysnewright)

    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)

    if maxindex != 0
        if maxvalue_type == 2 || maxvalue_type == 3
            Xarraysnew[maxindex] = constructoneXarray(Xpnew[maxindex],length(Xarraysnew[maxindex]),L)
        elseif maxvalue_type == 1
            Xarraysnew[maxindex_minus] = constructoneXarray(Xpnew[maxindex_minus],length(Xarraysnew[maxindex_minus]),L)
        else
            println("error in liquid merging, maxvalue_type is not 1,2,3")
        end
    end

    Xarraysnew,θarraysnew
end


function getinsertindex(Xp,Xvapornew_center,L,closedornot)
if closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center && Xp[index][2] >= Xp[index][1]
            return index
        end

        if (mod(Xvapornew_center - Xp[index][1],L) < mod(Xp[index][2] - Xp[index][1],L)) && Xp[index][2] <= Xp[index][1]
            return index
        end
    end
end


if !closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center
            return index
        end
    end
end

    return NaN
end

function getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)

    Xpnew = deepcopy(Xp)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)

    insertXp1=mod.((Xp[index][1]+Lliquid_adjust/2,Xvapornew[1]),L)
    insertXp2=mod.((Xvapornew[2],Xp[index][2]-Lliquid_adjust/2),L)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end

function getnewM(M,index,Minsert,closedornot)

    Mnew = deepcopy(M)

    insert!(Mnew,index+1,Minsert)

    return Mnew
end


function getsuperheat(Xstation,sys::UniformPhpSystem)
    @unpack PtoT = sys.propconvert
    Δθ = sys.mapping.θ_interp_walltoliquid(Xstation) - PtoT(sys.mapping.P_interp_liquidtowall(Xstation))
    return Δθ
end

function is_there_space_for_new_bubble(p::UniformPhpSystem,i)
    suitable_flag = true
    L_newbubble = p.wall.L_newbubble

    index = getinsertindex(_ξp(p),p.wall.Xstations[i],p.tube.L,p.tube.closedornot)

    L_liquid_left =  mod(p.wall.Xstations[i] - _ξp(p)[index][1],p.tube.L)
    L_liquid_right = mod(_ξp(p)[index][2] - p.wall.Xstations[i],p.tube.L)

    if (2*L_newbubble> L_liquid_left) || (2*L_newbubble > L_liquid_right)
        suitable_flag = false
    end

    return suitable_flag
end

function get_vapor_energy(sys0::UniformPhpSystem)
    volume_vapor = get_vapor_volume(sys0)
    PtoT = sys0.propconvert.PtoT
    return sum(PtoT.(sys0.vapor.P) .* sys0.liquid.Cpₗ .* volume_vapor)
end

function get_liquid_energy(sys0::UniformPhpSystem)
    Ac = sys0.tube.Ac
    ρ = sys0.liquid.ρₗ
    Cp = sys0.liquid.Cpₗ
    L = sys0.tube.L
    sys0.liquid.ξarrays;
    ds_liquid = [mod(Xarray[end]-Xarray[1],L)/length(Xarray) for Xarray in sys0.liquid.ξarrays];
    liquid_energy = sum.(sys0.liquid.θarrays).*ds_liquid .* (Ac .* ρ .* Cp)
    sum(liquid_energy)
end
