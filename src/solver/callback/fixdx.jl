export fixdx_affect!,fixdx_condition

# As liquid slugs grow and shrink, the internal conduction points
# (defined by ξarrays and θarrays) can get too bunched up or spread out.
# This callback watches out for those scenarios
# and rebuilds the internal conduction points if needed.

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

function getReconstructFlags(δξmin::Real,δξmax::Real,p::UniformPhpSystem)

    # only for closed loop tube
    numofliquidslug = length(p.liquid.Q_B)
    L = p.tube.L
    ReconstructFlags = Array{Bool,1}(undef, numofliquidslug)

    ξarrays = p.liquid.ξarrays

    for i in 1:numofliquidslug
        # Mark the internal conduction points of the i-th slug for reconstruction if 
        # a) grid spacing is less than a minimum value
        # b) grid spacing is greater than a maximum value
        ReconstructFlags[i] = mod((ξarrays[i][2] - ξarrays[i][1]),L) < δξmin || mod((ξarrays[i][2] - ξarrays[i][1]),L) > δξmax
    end

    return ReconstructFlags
end

function fixdx_affect!(integrator)
    p = deepcopy(getcurrentsys_nowall!(integrator.u,integrator.p));

    # support for nonuniform later
    @assert p isa UniformPhpSystem

    L = p.tube.L

    dξ_init = p.tube.L / p.tube.N

    δξmax = 2*dξ_init
    δξmin = 0.5*dξ_init

    sys = deepcopy(getcurrentsys_nowall!(integrator.u,integrator.p));

    ReconstructFlags = getReconstructFlags(δξmin,δξmax,p)
    indexReconstructSite = sort(findall(x->x == true, ReconstructFlags),rev = true)

    θ_interp_liquidtowall = p.mapping.θ_interp_liquidtowall

    ξ_L = p.liquid.ζ_L ./ p.tube.Ac
    ξ_R = p.liquid.ζ_R ./ p.tube.Ac

    Lliquid = ξ_to_liquid_length(ξ_L, ξ_R, sys.tube.L)
    Nliquid = round.(Int, Lliquid./dξ_init) .+ 1;
    # TODO: change after regression test
    # Nliquid = max.(2, round.(Int, Lliquid ./ dξ_init) .+ 1) # at least two internal conduction pts

    for i in indexReconstructSite
        p.liquid.ξarrays[i] = constructoneXarray((p.liquid.ξarrays[i][1],p.liquid.ξarrays[i][end]),Nliquid[i],L)
        p.liquid.θarrays[i] = θ_interp_liquidtowall(p.liquid.ξarrays[i])
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
