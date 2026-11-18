export vaporMergingAffect!,vaporMergingCondition,vaporMerging

function vaporMergingAffect!(integrator)

    p = deepcopy(getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p)));
    δv = DEFAULT_VAPOR_MERGE_FRAC*p.wall.L_newbubble


    merge_flags = getVaporMergeFlags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    for i in indexmergingsite
        p = vaporMerging(p,i)
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

function vaporMergingCondition(u,t,integrator)     # only for closed loop tube

    sys = getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p));
    δv = DEFAULT_VAPOR_MERGE_FRAC*sys.wall.L_newbubble
    merge_flags = getVaporMergeFlags(δv,sys)

    return sum(merge_flags) != 0
end

function vaporMerging(p::UniformPhpSystem,i)
    Xp = _ξp(p)

    Ac = p.tube.Ac
    ρₗ = p.liquid.ρₗ
    Mvapor_old = getMvapor(p);
    Mfilm_old = sum(getMfilm(p));
    Mliquid_old = getMliquid(p);

    Mold = sum(Mvapor_old + Mfilm_old + Mliquid_old)
    Lliquidslug_old = ξ_to_liquid_length(first.(Xp), last.(Xp), p.tube.L)
    Lvaporplug =  ξ_to_vapor_length(first.(Xp), last.(Xp),p.tube.L,p.tube.closedornot)
    left_index_after = [length(Lvaporplug)-1;1:(length(Lvaporplug)-1)]
    right_index_after = [1:length(Lvaporplug)-1;1]
    right_index = [2:length(Lvaporplug);1]


    systemp = deepcopy(p)
    Xpnew = _ξp(systemp)
    splice!(Xpnew,i)
    _set_ξp!(systemp, Xpnew)
    splice!(systemp.liquid.Q_B,i)
    splice!(systemp.liquid.ξarrays,i)
    splice!(systemp.liquid.θarrays,i)


    if i != length(Xp)
        splice!(systemp.vapor.Af_L,right_index[i])
        splice!(systemp.vapor.Af_R,i)
        splice!(systemp.vapor.Lfilm_start,right_index[i])
        splice!(systemp.vapor.Lfilm_end,i)
        splice!(systemp.vapor.P,i:i+1,[(systemp.vapor.P[i]+systemp.vapor.P[right_index[i]])/2])

        else
            splice!(systemp.vapor.Af_R,i)
            splice!(systemp.vapor.Lfilm_end,i)

            splice!(systemp.vapor.Af_L,length(systemp.vapor.Af_L));
            splice!(systemp.vapor.Af_L,1);
            insert!(systemp.vapor.Af_L,1,p.vapor.Af_L[end])

            splice!(systemp.vapor.Lfilm_start,length(systemp.vapor.Lfilm_start));
            splice!(systemp.vapor.Lfilm_start,1);
            insert!(systemp.vapor.Lfilm_start,1,p.vapor.Lfilm_start[end])

            splice!(systemp.vapor.P,length(systemp.vapor.P));
            splice!(systemp.vapor.P,1);
            insert!(systemp.vapor.P,1,(p.vapor.P[i]+p.vapor.P[right_index[i]])/2)
        end
    Mvapor_temp = getMvapor(systemp);
    Mfilm_temp = sum(getMfilm(systemp));
    Mliquid_temp = getMliquid(systemp);

    Mtemp = sum(Mvapor_temp + Mfilm_temp + Mliquid_temp)

    Mdiff = Mold - Mtemp

    @unpack PtoD = p.propconvert
    A = (Mdiff / Ac / (ρₗ - PtoD(systemp.vapor.P[right_index_after[i]])))
    B = (Lliquidslug_old[i] + p.vapor.Lfilm_end[i] + p.vapor.Lfilm_start[right_index[i]])

    L_lengthen =  A < B ? A : B

    Xp_new = _ξp(systemp)

    Xp_new[right_index_after[i]] = (Xp_new[right_index_after[i]][1] - L_lengthen/2,Xp_new[right_index_after[i]][2])
    Xp_new[left_index_after[i]] = (Xp_new[left_index_after[i]][1],Xp_new[left_index_after[i]][2] + L_lengthen/2)

    _set_ξp!(systemp, Xp_new)
    systemp.liquid.ξarrays = updateXarrays(Xp_new, systemp.liquid.θarrays, systemp.tube.L)

    return deepcopy(systemp)
end

function getVaporMergeFlags(δv::Real,sys::UniformPhpSystem)
    numofliquidslug = length(_ξp(sys))
    numofmergingsite = numofliquidslug
    merge_flags = Array{Bool,1}(undef, numofmergingsite)
    Xp = _ξp(sys)
    for i in 1:numofmergingsite
        merge_flags[i] = ((Xp[i][2] - Xp[i][1]) < δv && (Xp[i][2] - Xp[i][1]) >= 0) || ((Xp[i][2] - Xp[i][1] + sys.tube.L) < δv && (Xp[i][2] - Xp[i][1]) < 0) ? true : false
    end

    return merge_flags
end
