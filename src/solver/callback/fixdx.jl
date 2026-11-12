export fixdx_affect!,fixdx_condition
# boiling_condition,
function fixdx_condition(u,t,integrator)

        p = getcurrentsys_nowall!(integrator.u,integrator.p);

        dξ_init = p.tube.L / p.tube.N

        δξmax = 2*dξ_init
        δξmin = 0.5*dξ_init

        # sys = deepcopy(getcurrentsys_nowall!(integrator.u,integrator.p));

        ReconstructFlags = getReconstructFlags(δξmin,δξmax,p)

        return sum(ReconstructFlags) != 0
        # return true
end

function getReconstructFlags(δξmin,δξmax,p)

    # only for closed loop tube
    numofliquidslug = length(p.liquid.Xp)
    L = p.tube.L
    ReconstructFlags = Array{Bool,1}(undef, numofliquidslug)

    Xarrays = p.liquid.Xarrays

    for i in 1:numofliquidslug
     # merging bubble length threshold
        ReconstructFlags[i] = mod((Xarrays[i][2] - Xarrays[i][1]),L) < δξmin || mod((Xarrays[i][2] - Xarrays[i][1]),L) > δξmax ? true : false

    end

    return ReconstructFlags
end


function fixdx_affect!(integrator)

    p = deepcopy(getcurrentsys_nowall!(integrator.u,integrator.p));
    L = p.tube.L

    dξ_init = p.tube.L / p.tube.N

    δξmax = 2*dξ_init
    δξmin = 0.5*dξ_init

    sys = deepcopy(getcurrentsys_nowall!(integrator.u,integrator.p));


    ReconstructFlags = getReconstructFlags(δξmin,δξmax,p)
    indexReconstructSite = sort(findall(x->x == true, ReconstructFlags),rev = true)

    θ_interp_liquidtowall = p.mapping.θ_interp_liquidtowall

    Lliquid = ξ_to_liquid_length(first.(p.liquid.Xp), last.(p.liquid.Xp), sys.tube.L)
    Nliquid = round.(Int, Lliquid./dξ_init) .+ 1;

    for i in indexReconstructSite
        # println("reconstruct dx! in", i ," at ",integrator.t)

        p.liquid.Xarrays[i] = constructoneXarray((p.liquid.Xarrays[i][1],p.liquid.Xarrays[i][end]),Nliquid[i],L)
        p.liquid.θarrays[i] = θ_interp_liquidtowall(p.liquid.Xarrays[i])
    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)

    
    Ac = p.tube.Ac
    d = p.tube.Dh
    Af_L = p.vapor.Af_L
    Af_R = p.vapor.Af_R
    Lfilm_start = p.vapor.Lfilm_start
    Lfilm_end = p.vapor.Lfilm_end

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* Af_L - Lfilm_end .* Af_R
    @unpack PtoD = p.propconvert
    M = PtoD.(p.vapor.P) .* volume_vapor

    unew=[ζMAL_to_u(p.liquid.ζ_L,p.liquid.ζ_R,p.liquid.Q_B,M,Af_L,Af_R,Lfilm_start,Lfilm_end); liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end
