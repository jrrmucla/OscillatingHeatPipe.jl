export slugbc_affect!,slugbc_condition

function slugbc_condition(u,t,integrator)
    _reject_nonuniform(integrator.p)
    return true
end

function slugbc_affect!(integrator)
    sys = getcurrentsys_nowall!(integrator.u,_reject_nonuniform(integrator.p));
    unew = deepcopy(integrator.u)

    θarrays = sys.liquid.θarrays

    adjacentT = getadjacentT.([sys],1:length(θarrays))

    for i in eachindex(θarrays)
        θarrays[i][1] = adjacentT[i][1]
        θarrays[i][end] = adjacentT[i][end]
    end

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), unew)

    unew[index_dynamics_end:end] = liquidθtovec(θarrays)

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end


function getadjacentT(p::UniformPhpSystem,i::Int64)
    @unpack PtoT = p.propconvert
    Tfirst = PtoT(p.vapor.P[i])

    if p.tube.closedornot == true
        Tlast = i >= length(p.vapor.P) ? PtoT(p.vapor.P[1]) : PtoT(p.vapor.P[i+1])
    else
        Tlast = PtoT(p.vapor.P[i+1])
    end

    (Tfirst,Tlast)
end
