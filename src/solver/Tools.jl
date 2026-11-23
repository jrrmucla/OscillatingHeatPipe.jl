export getgvec,getgh,delta, # get actrual heightg of the tube
ξ_to_vapor_length,ξ_to_liquid_length,ζ_to_vapor_volume,ζ_to_liquid_volume,getXpvapor, # transfer Xp to the length of vapors, length of liquids, and Xp for vapor.
# ifamong,constructXarrays,
Hfilm,getδarea,getδFromδarea,getMvapor,getMfilm,getMliquid,getMtotal,getchargeratio,get_vapor_volume,
getCa,f_churchill,Catoδ,RntoΔT

# typing of NamedTuple as an alternative to AbstractPHP
# is to allow for unit tests with minimal PHP-like tuples

function getgvec(g0::T,g_angle::T=3/2*π) where {T<:Real}
    g = g0*[cos(g_angle),sin(g_angle)]
end

function getgh(g::Vector{T},x::Vector{T},y::Vector{T}) where {T<:Real}
    xy = [x';y']
    vec(sum(-g .* xy,dims=1));
end

delta(n) = n != 0

"""
    ξ_to_vapor_length(ξ_L, ξ_R, L, closedornot)

Return vapor plug lengths from interface ξ-locations.
- closed tubes: plug `i` lies between slugs `i-1` and `i`
- open tubes: plug `i` lies between slugs `i` and `i+1`
"""
function ξ_to_vapor_length(ξ_L::AbstractVector{<:Real}, ξ_R::AbstractVector{<:Real}, L::Real, closedornot::Bool)
    length(ξ_L) == length(ξ_R) ||
        throw(DimensionMismatch("Left and right interface counts must match"))
    if closedornot
        isempty(ξ_L) && throw(ArgumentError("A closed tube requires at least one liquid slug"))
        return mod.(ξ_L - circshift(ξ_R, 1), L)
    end
    length(ξ_L) > 1 || throw(ArgumentError("An open tube requires at least two liquid slugs"))
    return ξ_L[2:end] - ξ_R[1:end-1]
end
"""
    ξ_to_liquid_length(ξ_L, ξ_R, L)

Return liquid slug lengths from interface ξ-locations.
"""
function ξ_to_liquid_length(ξ_L::Real, ξ_R::Real, L::Real)
    return mod(ξ_R - ξ_L, L)
end

function ξ_to_liquid_length(ξ_L::AbstractVector{<:Real}, ξ_R::AbstractVector{<:Real}, L::Real)
    return mod.(ξ_R - ξ_L, L)
end

"""
    ζ_to_liquid_volume(ζ_L, ζ_R, V)

Return liquid slug volumes from interface ζ-locations.
"""
ζ_to_liquid_volume(ζ_L::Real, ζ_R::Real, V::Real) = ξ_to_liquid_length(ζ_L, ζ_R, V)
ζ_to_liquid_volume(ζ_L::AbstractVector{<:Real}, ζ_R::AbstractVector{<:Real}, V::Real) =
    ξ_to_liquid_length(ζ_L, ζ_R, V)

"""
    ζ_to_vapor_volume(ζ_L, ζ_R, V, closedornot)

Return vapor plug volumes from interface ζ-locations.
- closed tubes: plug `i` lies between slugs `i-1` and `i`
- open tubes: plug `i` lies between slugs `i` and `i+1`
These volumes include the liquid films.
"""
ζ_to_vapor_volume(ζ_L::AbstractVector{<:Real}, ζ_R::AbstractVector{<:Real}, V::Real, closedornot::Bool) =
    ξ_to_vapor_length(ζ_L, ζ_R, V, closedornot)

"""
    getXpvapor(Xp::Vector{Tuple},closedornot) -> Vector{Tuple}

Given the vector of liquid slug coordinates `Xp`, return the vector of vapor plug coordinates. 
"""
function getXpvapor(Xp::Vector{Tuple{Float64, Float64}},closedornot::Bool)

    # this is not periodic -- it has ends.  
    if closedornot == false 
        Xpvapor = map(tuple,zeros(length(Xp)-1),zeros(length(Xp)-1))
        maxindex = length(Xp) - 1
        if maxindex > 1
            for i = 1:maxindex
                Xpvapor[i] = (Xp[i][end], Xp[i+1][1])
            end
        else
            return error("need more than one liquid slugs")
        end
    end

    # this is a periodic tube, so there are as many vapor plugs and liquid slugs
    if closedornot == true  
        Xpvapor = map(tuple,zeros(length(Xp)),zeros(length(Xp)))
        Xpvapor[1]=(Xp[end][end],Xp[1][1])
        for i = 2:(length(Xp))
            Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
        end
    end

    return Xpvapor
end


"""
    ifamongone(value::Real,range::Tuple)

Determine if `value` lies between the pair of values in the tuple `range`.
This accounts for cases in closed tubes in which `range` crosses the branch cut,
so that the start point in `range` has a larger value than the end point.
"""
function ifamongone(value::Float64, range::Tuple{Float64,Float64})
    # Note: this function does not make use of the tube length, so it might have errors.
    # It should also be possible to do it faster.
    return ((value >= range[1]) && (value <= range[end])) || ((value <= range[end]) && (range[1] >= range[end])) || ((value >= range[1]) && (range[1] >= range[end])) ? true : false
end

"""
    ifamong(value::Real,X::Vector{Tuple})

Determine if `value` lies between any of the tuples in the vector `X`.
"""
function ifamong(value::Float64, X::Vector{Tuple{Float64,Float64}})
    return Bool(sum(ifamongone.(value,X)) >= 1 ? true : false)
end

"""
    constructXarrays(X::Vector{Tuple},N::Int,θi::Float64,L::Float64) -> Vector{Vector}, Vector{Vector}

Initialize the 1D arc coordinate and temperature fields for every liquid slug. The
number of field points in each slug is never smaller than 2, and set equal to the fraction
of `N` based on the slug's length relative to `L`. 
The coordinate array is interpolated between the values in each `X` element (with 
the branch cut at `L` in a closed tube respected), and the temperature array
is set to `θi` uniformly. 
"""
function constructXarrays(X0::Vector{Tuple{Float64, Float64}},N::Int,θinitial::Float64,L::Float64)
    Xarrays=Array{Array{Float64, 1}, 1}(undef, length(X0))

    Lliquid = ξ_to_liquid_length(first.(X0), last.(X0), L)
    Nliquid = max.(ceil.(Int, N.*Lliquid./L),2)

    for i = 1:length(Xarrays)
        if X0[i][1] < X0[i][2]
            Xarrays[i] = range(X0[i][1], X0[i][2], length=Nliquid[i])
        else
            Xarrays[i] = mod.(range(X0[i][1], X0[i][2]+L, length=Nliquid[i]), L)
            # Xarrays[i] = mod.(Xarrays[i], L)
        end
    end

    θarrays = deepcopy(Xarrays)
    for i in eachindex(θarrays)
        θarrays[i][:] .= θinitial
    end

    return(Xarrays,θarrays)
end

"""
    constructoneXarray(X::Tuple,N::Int,L::Float64) -> Vector

Initialize the 1D arc coordinate and temperature fields for one liquid slug. The
number of field points in each slug is set equal to `N`. 
The coordinate array is interpolated between the values in `X` (with 
the branch cut at `L` in a closed tube respected), and the temperature array
is set to `θi` uniformly. 
"""
function constructoneXarray(X0::Tuple{Float64,Float64},Nliquid,L)
    Xarray=Array{Float64, 1}

        if X0[1] < X0[2]
            Xarray = Array(range(X0[1], X0[2], length=Nliquid))
        else
            Xarray = Array(mod.(range(X0[1], X0[2]+L, length=Nliquid),L))
        end

        Xarray[1] = X0[1]
        Xarray[end] = X0[2]

    return Xarray
end


"""
    initialize X and θ field for wall, return Array{Float64, 1} and Array{Float64, 1}

    X0       :: Array{Tuple{Float64,Float64},1}
    N        :: Int, the number of cells in the wall (ΔX for liquid equals ΔX for the wall)
    θinitial :: value
    L        :: tube length
"""

"""
    constructwallXθarray(X::Vector,θi::Float64,curv::Vector) -> Vector, Vector, Vector

Return wall vectors for arc length coordinate, temperature, and curvature, given
the vector `X` for coordinates, a temperature value `θi`, and a vector of curvatures
`curv` 
"""
function constructwallXθarray(line::Vector{Float64},θinitial::Float64,curv::Vector{Float64})
    Xwallarray = deepcopy(line)
    θwallarray = Xwallarray .* 0 .+ θinitial
    curvwallarray = deepcopy(curv)

    return(Xwallarray,θwallarray,curvwallarray)
end



"""
    Hfilm(δfilm,sys::UniformPhpSystem) -> Float64

Return the heat transfer coefficient for the film with thickness `δfilm`
"""
Hfilm(δfilm,sys::Union{AbstractPHP, NamedTuple}) = Hfilm(δfilm,sys.vapor)

function Hfilm(δfilm,vapor::Vapor)
    @unpack δmin,k,Hᵥ  = vapor
    δthreshold = FILM_THRESHOLD
    δmax = FILM_MAX_THICKNESS
    kₗ = k

    if (δfilm > δthreshold) && (δfilm <= δmax)
        return kₗ/δfilm
    elseif (δfilm > δmax) && (δfilm <= 2δmax)
        return  kₗ/δmax - (δfilm-δmax)*(kₗ/δmax^2) + 1e-6
    elseif δfilm > δmin && δfilm <= δthreshold
        return  Hᵥ + (δfilm-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin) + 1e-6
    else
        # return Hᵥ  + 1e-6
        return 0.0
    end
end



"""
    getδarea(Ac,d,δ)

Given cross-sectional tube area `Ac`, diameter `d`, and film thickness `δ`,
return `Ac * (1 - (1 - 2δ/d)^2)`. This is the exact annular area for a
circular channel, with `0 ≤ δ ≤ d/2`. It also describes a uniform wall film
in a square channel when `d` is its side length; arbitrary cross sections
require their own geometry relation.
"""
getδarea(Ac::Real,Dh::Real,δ::Union{Real,AbstractArray{<:Real}}) = Ac .* (2δ ./ Dh) .* (2 .- 2δ ./ Dh)

"""
    getδFromδarea(Ac,d,δarea)

Given cross-sectional tube area `Ac`, diameter `d`, and film cross-sectional
area `δarea`, return the film thickness `d/2 * (1 - sqrt(1 - δarea/Ac))`.
For circular channels this exactly inverts the annular area; `0 ≤ δarea ≤ Ac`.
"""
getδFromδarea(Ac,d,δarea) = (d / 2) * (δarea / Ac) / (1 + sqrt(1 - δarea / Ac))
# Rationalized inverse of 1 - sqrt(1 - Afilm/Ac), avoiding cancellation
# for very thin films. Invalid areas above Ac retain a DomainError.

"""
Return the masses of all of the vapor,liquid, and film regions
"""
function getMtotal(sys::Union{AbstractPHP, NamedTuple})
    return sum(getMvapor(sys)) + sum(getMliquid(sys)) + sum(sum(getMfilm(sys)))
end

"""
    getchargeratio(sys::UniformPhpSystem) -> Scalar

Return the chargeratio considering all of the vapor,liquid, and film regions

TODO: this isn't consistent with initialization scheme
and it looks like initialization is in the wrong
"""
function getchargeratio(sys::Union{AbstractPHP, NamedTuple})
    Mtotal = getMtotal(sys)
    @unpack ρₗ = sys.liquid # get the liquid density
    Vtotal = _tube_volume(sys.tube)
    Mliquidfull = ρₗ * Vtotal # mass of the tube if it is full of liquid
    return Mtotal / Mliquidfull
end
 
"""
    getMvapor(sys::UniformPhpSystem) -> Vector

Return the masses of all of the vapor regions
"""
function getMvapor(sys::Union{AbstractPHP, NamedTuple})

    @unpack propconvert, vapor = sys
    @unpack PtoD = propconvert
    @unpack P = vapor

    ρᵥ = PtoD.(P)

    return ρᵥ .* get_vapor_volume(sys)
    
end


"""
    get_vapor_volume(sys)

Return the volumes of the vapor regions, given the tube system `sys`.
"""
function get_vapor_volume(sys::Union{AbstractPHP, NamedTuple})
    @unpack tube, vapor, liquid = sys
    @unpack Af_L, Af_R, Lfilm_start, Lfilm_end = vapor
    @unpack ζ_L, ζ_R = liquid
    plug_volumes = ζ_to_vapor_volume(ζ_L, ζ_R, _tube_volume(tube), tube.closedornot)
    return plug_volumes .- Af_L .* Lfilm_start .- Af_R .* Lfilm_end
end

"""
Return the masses of all of the liquid films in the system.
"""
function getMfilm(sys::Union{AbstractPHP, NamedTuple})
    @unpack tube, vapor, liquid = sys
    @unpack Af_L, Af_R, Lfilm_start, Lfilm_end = vapor
    @unpack ρₗ = liquid

    return _getM.(ρₗ,Af_L,Lfilm_start), _getM.(ρₗ,Af_R,Lfilm_end)
end

"""
Return the masses of all of the liquid slugs in the system.
"""
function getMliquid(sys::Union{AbstractPHP, NamedTuple})
    @unpack tube, liquid = sys
    @unpack ρₗ, ζ_L, ζ_R = liquid
    volumes = tube.closedornot ? ζ_to_liquid_volume(ζ_L, ζ_R, _tube_volume(tube)) : ζ_R - ζ_L
    return ρₗ .* volumes
end

_tube_region_volume(tube::UniformTube, start, region_length) = tube.Ac * region_length

function _tube_region_volume(tube::VariableTube, start, region_length)
    tube_length = _tube_length(tube)
    stop = start + region_length
    if stop > tube_length
        return _tube_volume_between(tube, start, tube_length) +
               _tube_volume_between(tube, 0.0, stop - tube_length)
    end
    return _tube_volume_between(tube, start, stop)
end

function _tube_volume_between(tube::VariableTube, start, stop)
    knots = tube.ξ_geometric_knots
    interior = knots[(knots .> start) .& (knots .< stop)]
    positions = [start; interior; stop]
    area = construct_geometric_interpolations(tube)[3].(positions)
    # trapezoidal integration is exact for piecewise-linear tube area
    return sum(diff(positions) .* (area[1:end-1] .+ area[2:end]) ./ 2)
end

_getM(ρₗ,A,L) = ρₗ*A*L


"""
    getCa(μ,σ,velocity)

Return the capillary numbers for all liquid interfaces
"""
function getCa(μ,σ,velocity)
    Ca = abs.(μ.*velocity./σ)
end

"""
    f_churchill(Re,ϵ)

Return the Darcy friction factor, based on the Churchill correlation.
"""
function f_churchill(Re,eps_rel=0.001)
    Θ1 = (-2.457*log((7/Re)^(0.9)  +  0.27 * eps_rel))^16
    Θ2 = (37530/Re)^16
    f=8*((8/Re)^12+(1/(Θ1+Θ2)^1.5))^(1/12)
    
    f
end

"""
    Catoδ(d,Ca[;adjust_factor=1,δmin=2e-6,δmax=1e-4])

Given capillary number `Ca` and `d`, return the thickness.
This is a generalization of the Aussillous and Quere formula.
It uses an adjustment factor that can be tuned empirically.
"""
function Catoδ(d::Union{AbstractFloat,AbstractVector},Ca;adjust_factor=1,δmin=2e-6,δmax=1e-4)
    δ = 1.34 .* Ca .^ (2/3) ./ (1 .+ 3.35 .* Ca .^ (2/3)) .* d ./ 2 .* adjust_factor
    return clamp.(δ, δmin, δmax)
end

"""
    RntoΔT(Rn,Tavg,fluid_type,d,TtoP) -> Vector

Given the nucleation radius `Rn`, the fluid type `fluid_type` and reference
temperature `Tavg`, the channel diameter `d`, and the temperature-to-pressure
relation `TtoP`, return the superheat threshold temperature for boiling.
"""
function RntoΔT(Rn,Tavg,fluid_type,d,TtoP)
    p_fluid = SaturationFluidProperty(fluid_type,Tavg);

    Rkg = p_fluid.R/p_fluid.M
    Rin = d/2
    P = TtoP(Tavg)

    y = Rkg .* Tavg ./ (p_fluid.hᵥ-p_fluid.hₗ) .* log.(1 .+ 2 .* p_fluid.σ ./ P .* (1 ./ Rn .- 1/(2Rin)))
    ΔTavg = Tavg .* (1 ./ (1 .- y) .- 1)
end

"""
    systoM(sys::AbstractPHP)

Given tube system `sys`, return the masses of the vapor regions.
"""
function systoM(sys0::AbstractPHP)
    return getMvapor(sys0)
end
