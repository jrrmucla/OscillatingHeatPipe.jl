export duliquidθtovec,duwallθtovec,liquidθtovec,wallθtovec,newstate,getcurrentsys_nowall!,getcurrentsys_nowall

#=
Some notes on structures of liquid and vapor positions:
- The Xp is coupled by every liquid slug. For instance, if there is one liquid slug. Xp is a one-element tuple (Xp[1][1], Xp[1][2]).
    But sometimes we need Xp to be coupled by every vapor plug. For one liquid slug, we have two vapor plugs.
    So by adding 0 and L at the beginning and the end,
    we construct a two-element tuple ((0.0,Xp[1][1]) and ((Xp[1][2],L). Generally, for every N-element Xp, we construct an N+1 element Xpvapor
        Xp    ::   the locations of all interfaces, each element means a liquid slug.
        L     ::   the length of the 1D tube
- Xp contains the liquid slug begin/end positions (in arclength). It is arranged as a
  vector of tuples, where the first tuple element is the beginning and second element is
  the end
- If the tube is periodic, then `tube.closedornot = true` and there are as many vapor plugs as liquid slugs.
  Vapor plug i is bounded by liquid slugs i-1 and i. Plug 1 is bounded by N and 1.
- If the tube is not periodic, then there is one less vapor plug than liquid slugs. Vapor plug i is bounded
  by slugs i and i+1.
=#

"""
    X_to_u(Xp,dXdt) -> Vector

Transform vectors `Xp`, `dXdt` with the coordinates and velocities of the liquid slug interfaces into
a single state vector 
"""
function X_to_u(Xp::Array{Tuple{Float64,Float64},1},dXdt::Array{Tuple{Float64,Float64},1})
    Np = length(Xp)
    u = zeros(4*Np)
    
    for i = 1:Np
        # input Xp
        u[2*i-1] = Xp[i][1]
        u[2*i] = Xp[i][end]
        # input dXdt
        u[2*Np + 2*i-1] = dXdt[i][1]
        u[2*Np + 2*i] = dXdt[i][end]
    end
    
    return u
end


"""
    XMAL_to_u(Xp,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end) -> Vector

Assemble vectors `Xp`, `dXdt` of the coordinates and velocities of the slug interfaces, and
mass `M` of the vapor regions, `Afilm_start` and `Afilm_end` of cross-sectional film areas (m²), and `Lfilm_start`
and `Lfilm_end` of film lengths, to form state vector.

The blocks are `[positions; velocities; masses; Afilm_start; Afilm_end; Lfilm_start; Lfilm_end]`.
Film areas are in m²; derive thickness with `getδFromδarea(Ac, d, A)`.
"""
function XMAL_to_u(Xp,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end)
    return ([X_to_u(Xp,dXdt);M;Afilm_start;Afilm_end;Lfilm_start;Lfilm_end])
end

"""
    u_to_XMAL(u) -> Vector{Tuple}, Vector{Tuple}, Vector, Vector, Vector, Vector, Vector

Disassemble the state vector `u` into vectors `Xp`, `dXdt` of the coordinates and velocities of the slug interfaces, and
mass `M` of the vapor regions, `Afilm_start` and `Afilm_end` of cross-sectional film areas (m²), and `Lfilm_start`
and `Lfilm_end` of film lengths.
"""
function u_to_XMAL(u::Array{Float64,1})

    maxindex = div(length(u),9)
    maxliquidindex = 0
    if mod(length(u),9) == 4
        maxliquidindex = maxindex + 1
    elseif mod(length(u),9) == 0
        maxliquidindex = maxindex
    else
        println("new function error")
        return error("new function error")
    end

    Xp   = map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    dXdt = map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    
    for i = 1:maxliquidindex
        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])
        # input dXdt
        dXdt[i] = (u[2*maxliquidindex + 2*i-1],u[2*maxliquidindex + 2*i])
    end

    M =             u[4*maxliquidindex+1:           4*maxliquidindex+maxindex]
    Afilm_start =        u[4*maxliquidindex+maxindex+1:  4*maxliquidindex+2maxindex]
    Afilm_end =          u[4*maxliquidindex+2maxindex+1: 4*maxliquidindex+3maxindex]
    Lfilm_start =   u[4*maxliquidindex+3maxindex+1: 4*maxliquidindex+4maxindex]
    Lfilm_end =     u[4*maxliquidindex+4maxindex+1: 4*maxliquidindex+5maxindex]
    return Xp,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end
end

"""
    A bunch of functions to transfer θ to state vector rate du
"""

function duliquidθtovec(duθarrays::Vector{Vector{Float64}})
    return vcat(map(oneduliquidθtovec, duθarrays)...)
end

function oneduliquidθtovec(duθarray::Vector{Float64})
    return [0.0; duθarray]
end

"""
    A bunch of functions to transfer θ to state vector u
"""

function liquidθtovec(θarrays::Vector{Vector{Float64}})
    return vcat(map(wallθtovec, θarrays)...)
end

function wallθtovec(θwall::Vector{Float64})
    return [SEPERATION_VAR; θwall]
end

function newstate(sys0::UniformPhpSystem)
    @unpack Ac,d,L,closedornot = sys0.tube
    @unpack Xp,dXdt = sys0.liquid
    @unpack Afilm_start,Afilm_end,Lfilm_start,Lfilm_end,P = sys0.vapor
    @unpack PtoD = sys0.propconvert

    M = systoM(sys0::UniformPhpSystem)

    u=[XMAL_to_u(Xp,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
end

"""
    Given a new u and an old system, return a new system `sysnew`.
    Does not update wall temperature. 
    Updates interpolations based on the current wall temperature in sys0.

    So this function is basically updating the system based on the most current u and wall temperature in sys0.
"""
function getcurrentsys_nowall!(u,sys0)

    indexes = Int64[]
    θliquidrec = Array[]

    for i in eachindex(u)
        if abs(u[i]+1e10) <= 10^(-1)
            push!(indexes,i)
        end
    end

    Xp,dXdt,M,Afilm_start,Afilm_end,Lfilm_start,Lfilm_end = u_to_XMAL(u[1:indexes[1]-1])

    if sys0.tube.closedornot == true 
        modX!(Xp,sys0.tube.L)
    end
    for i = 1:length(indexes)-1
    push!(θliquidrec, u[indexes[i]+1:indexes[i+1]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L,sys0.tube.closedornot)

    Ac = sys0.tube.Ac

    d = sys0.tube.d

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* Afilm_start - Lfilm_end .* Afilm_end
    ρ = M ./ volume_vapor

    # println(Lvaporplug)
    # println(Lfilm_start)
    # println(Lfilm_end)

    @unpack DtoP = sys0.propconvert
    P = DtoP.(ρ)
  
    sysnew = sys0

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.Xarrays = updateXarrays(Xp,sysnew.liquid.θarrays,sysnew.tube.L)

    sysnew.vapor.P = P
    sysnew.vapor.Afilm_start = Afilm_start
    sysnew.vapor.Afilm_end = Afilm_end
    sysnew.vapor.Lfilm_start = Lfilm_start
    sysnew.vapor.Lfilm_end = Lfilm_end

    sysnew.mapping = refresh_interpolations(sysnew)

    return sysnew
end

