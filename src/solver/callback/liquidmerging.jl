export merging_affect!,merging_condition,nucleateboiling,merging

function merging_affect!(integrator)

    p = deepcopy(getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p)));
    δv = DEFAULT_LIQUID_MERGE_FRAC*p.wall.L_newbubble

    merge_flags = getmerge_flags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    for i in indexmergingsite
        p = merging(p,i)
    end



    Af_L = p.vapor.Af_L
    Af_R = p.vapor.Af_R
    Lfilm_start = p.vapor.Lfilm_start
    Lfilm_end = p.vapor.Lfilm_end

    volume_vapor = get_vapor_volume(p)
    @unpack PtoD = p.propconvert
    M = PtoD.(p.vapor.P) .* volume_vapor

    unew=[ζMAL_to_u(p.liquid.ζ_L,p.liquid.ζ_R,p.liquid.Q_B,M,Af_L,Af_R,Lfilm_start,Lfilm_end); liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end

function merging_condition(u,t,integrator)     # only for closed loop tube

    p = getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p));
    δv = DEFAULT_LIQUID_MERGE_FRAC*p.wall.L_newbubble

    sys = p

    merge_flags = getmerge_flags(δv,sys)

    return sum(merge_flags) != 0
end

function merging(sys::UniformPhpSystem,i)
    Xp = _ξp(sys)

    closedornot = sys.tube.closedornot
    Lliquidslug = ξ_to_liquid_length(first.(Xp), last.(Xp), sys.tube.L)
    Lvaporplug  = ξ_to_vapor_length(first.(Xp), last.(Xp),sys.tube.L,sys.tube.closedornot)

    numofliquidslug = length(Lliquidslug)
    numofvaporbubble = length(Lvaporplug)


    L = sys.tube.L
    Ac = sys.tube.Ac
    ρₗ = sys.liquid.ρₗ

    Mvapor = getMvapor(sys)
    Mfilm = getMfilm(sys);

    @unpack PtoD = sys.propconvert
    if closedornot == true
        left_index = i > 1 ? i-1 : length(Lvaporplug)
        right_index = i < length(Lvaporplug) ? i+1 : 1
        liquid_right_i = i #  operating liquid i
        liquid_left_i = left_index

        splice_liquid_i_1 = i > 1 ? left_index : numofliquidslug
        splice_liquid_i_2 = i > 1 ? left_index : 1
        insert_i_1 = i > 1 ? left_index : numofliquidslug-1

    else
        left_index = i-1
        right_index = i+1
        liquid_right_i = i+1 #  operating liquid i
        liquid_left_i = i

        splice_liquid_i_1 = i
        splice_liquid_i_2 = i
        insert_i_1 = i


    end



    if closedornot == false && i == 1
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(sys.vapor.P[right_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(sys.vapor.P[right_index])))
        Xpnewone = (Xp[liquid_left_i][1], Xp[liquid_right_i][end] - Lvaporplug[i] + Linsert)
        Q_B_new = 0.0
    elseif closedornot == false && i == numofvaporbubble
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(sys.vapor.P[left_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(sys.vapor.P[left_index])))
        Xpnewone = (Xp[liquid_left_i][1]+Lvaporplug[i] - Linsert, Xp[liquid_right_i][end])
        Q_B_new = 0.0
    else
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(sys.vapor.P[left_index]) .+ PtoD(sys.vapor.P[right_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(sys.vapor.P[left_index]) .+ PtoD(sys.vapor.P[right_index])))
        Xpnewone = (mod(Xp[liquid_left_i][1]+Lvaporplug[i]/2 - Linsert/2,L), mod(Xp[liquid_right_i][end] - Lvaporplug[i]/2 + Linsert/2,L))
        Q_B_new = (sys.liquid.Q_B[liquid_left_i]*Lliquidslug[liquid_left_i] + sys.liquid.Q_B[liquid_right_i]*Lliquidslug[liquid_right_i])/(Lliquidslug[liquid_left_i]+Lliquidslug[liquid_right_i])
    end

    sys_temp = deepcopy(sys)

    Xp_new = _ξp(sys_temp)
    splice!(Xp_new,splice_liquid_i_1)
    splice!(Xp_new,splice_liquid_i_2)
    insert!(Xp_new,insert_i_1,Xpnewone)
    _set_ξp!(sys_temp, Xp_new)

    splice!(sys_temp.liquid.Q_B,splice_liquid_i_1)
    splice!(sys_temp.liquid.Q_B,splice_liquid_i_2)
    insert!(sys_temp.liquid.Q_B,insert_i_1,Q_B_new)

    splice!(sys_temp.vapor.Af_L,i)
    splice!(sys_temp.vapor.Af_R,i)
    splice!(sys_temp.vapor.Lfilm_start,i)
    splice!(sys_temp.vapor.Lfilm_end,i)
    splice!(sys_temp.vapor.P,i)
    Nliquids = length([sys_temp.liquid.ξarrays[liquid_left_i]; sys_temp.liquid.ξarrays[liquid_right_i]])
    Xarraysnewone = constructoneXarray(Xp_new[insert_i_1],Nliquids-1,sys.tube.L)
    splice!(sys_temp.liquid.ξarrays,splice_liquid_i_1);
    splice!(sys_temp.liquid.ξarrays,splice_liquid_i_2);
    insert!(sys_temp.liquid.ξarrays,insert_i_1,Xarraysnewone);
    θarraysnewone = [sys.liquid.θarrays[liquid_left_i][1:end-1]; (sys.liquid.θarrays[liquid_left_i][end-1]+sys.liquid.θarrays[liquid_right_i][2])/2 ;sys.liquid.θarrays[liquid_right_i][2:end]]
    splice!(sys_temp.liquid.θarrays,splice_liquid_i_1);
    splice!(sys_temp.liquid.θarrays,splice_liquid_i_2);
    insert!(sys_temp.liquid.θarrays,insert_i_1,θarraysnewone);

    return deepcopy(sys_temp)
end

function getmerge_flags(δv::Real,sys::UniformPhpSystem)

    Xpvapor = getXpvapor(_ξp(sys),sys.tube.closedornot)
    numofmergingsite = length(Xpvapor)
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    for i in 1:numofmergingsite
        merge_flags[i] = sys.tube.closedornot ? (mod(Xpvapor[i][2] - Xpvapor[i][1],sys.tube.L) < δv) : ((Xpvapor[i][2] - Xpvapor[i][1]) < δv)
    end

    return merge_flags
end
