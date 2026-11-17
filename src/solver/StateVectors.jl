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
    ζMAL_to_u(ζ_L, ζ_R, Q_B, M, Af_L, Af_R, Lfilm_start, Lfilm_end)

Assemble `[ζ_L; ζ_R; Q_B; M; Af_L; Af_R; Lfilm_start; Lfilm_end]`.
Interface coordinates are volumes (m³), bulk flow rates are m³/s, and film areas are m².
"""
function ζMAL_to_u(ζ_L::Array{Float64,1}, ζ_R::Array{Float64,1}, Q_B::Array{Float64,1},M,Af_L,Af_R,Lfilm_start,Lfilm_end)
    return ([ζ_L;ζ_R;Q_B;M;Af_L;Af_R;Lfilm_start;Lfilm_end])
end

"""
    u_to_ζMAL(u)

Return the eight blocks `(ζ_L, ζ_R, Q_B, M, Af_L, Af_R, Lfilm_start, Lfilm_end)`.
Closed tubes have equal liquid/vapor counts; open tubes have one extra liquid slug.
"""
function u_to_ζMAL(u::AbstractVector)
    # closed case: length = 8nL
    # open case: length = 8nL-5
    nV = div(length(u),8)

    if mod(length(u),8) == 3
        nL = nV + 1 # open ohp
    elseif mod(length(u),8) == 0
        nL = nV # closed ohp
    else
        return error("new function error")
    end

    ζ_L =   u[1:nL]
    ζ_R =   u[nL+1:2nL]
    Q_B =   u[2nL+1:3nL]
    M   =   u[3nL+1:3nL+nV]
    Af_L =          u[3nL+nV+1:3nL+2nV]
    Af_R =          u[3nL+2nV+1:3nL+3nV]
    Lfilm_start =   u[3nL+3nV+1:3nL+4nV]
    Lfilm_end =     u[3nL+4nV+1:3nL+5nV]
    return ζ_L,ζ_R,Q_B,M,Af_L,Af_R,Lfilm_start,Lfilm_end
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
    @unpack ζ_L, ζ_R, Q_B = sys0.liquid
    @unpack Af_L,Af_R,Lfilm_start,Lfilm_end,P = sys0.vapor
    @unpack PtoD = sys0.propconvert

    M = systoM(sys0::UniformPhpSystem)

    u=[ζMAL_to_u(ζ_L,ζ_R,Q_B,M,Af_L,Af_R,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
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

    ζ_L,ζ_R,Q_B,M,Af_L,Af_R,Lfilm_start,Lfilm_end = u_to_ζMAL(u[1:indexes[1]-1])

    if sys0.tube.closedornot == true 
        tube_volume = _tube_volume(sys0.tube)
        ζ_L .= mod.(ζ_L, tube_volume)
        ζ_R .= mod.(ζ_R, tube_volume)
    end
    for i = 1:length(indexes)-1
    push!(θliquidrec, u[indexes[i]+1:indexes[i+1]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    ξ_L = map_ζ_to_ξ(sys0.tube, ζ_L)
    ξ_R = map_ζ_to_ξ(sys0.tube, ζ_R)
    Xp = collect(zip(ξ_L, ξ_R))

    sysnew = sys0

    sysnew.liquid.ζ_L = ζ_L
    sysnew.liquid.ζ_R = ζ_R
    sysnew.liquid.Q_B = Q_B
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.ξarrays = updateXarrays(Xp,sysnew.liquid.θarrays,sysnew.tube.L)

    sysnew.vapor.Af_L = Af_L
    sysnew.vapor.Af_R = Af_R
    sysnew.vapor.Lfilm_start = Lfilm_start
    sysnew.vapor.Lfilm_end = Lfilm_end

    volume_vapor = get_vapor_volume(sysnew)
    ρ = M ./ volume_vapor
    @unpack DtoP = sysnew.propconvert
    sysnew.vapor.P = DtoP.(ρ)

    sysnew.mapping = refresh_interpolations(sysnew)

    return sysnew
end

