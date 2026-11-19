export PhpSystem_nomapping,UniformPhpSystem,Liquid,Vapor,Wall,PropConvert,Mapping,Cache

using Parameters

struct PropConvert
    fluid_type::String
    PtoT    :: AbstractInterpolation
    TtoP    :: AbstractInterpolation
    PtoD    :: AbstractInterpolation
    DtoP    :: AbstractInterpolation
    PtoHfg  :: AbstractInterpolation
end

function PropConvert(fluid_type::String);
    PtoT,TtoP,PtoD,DtoP,PtoHfg = createCoolPropinterpolation(fluid_type::String)
    PropConvert(fluid_type,PtoT,TtoP,PtoD,DtoP,PtoHfg);
end

"""
Liquid is a struct containing liquid properties at a ref temperature
    Hₗ::Float64                              heat transfer coefficient between the wall and the pure liquid slug
    ρ::Float64                              liquid density
    Cp::Float64                             liquid specific heat capacity
    α::Float64                              liquid heat diffusivity
    μₗ::Float64                              liquid dynamic viscosity
    σ::Float64                              liquid surface tension
    ζ_L::Array{Tuple{Float64,Float64},1}     location (in volume coordinates) of the left interface of the i-th slug
    ζ_R::Array{Tuple{Float64,Float64},1}     location (in volume coordinates) of the right interface of the i-th slug
    Q_B::Array{Tuple{Float64,Float64},1}     bulk volumetric flow rate of the i-th slug
    ξarrays::Array{Array{Float64,1},1}      finite difference location points (in length coordinates) within each liquid slug
    θarrays::Array{Array{Float64,1},1}      finite difference temperature points within each liquid slug
"""
mutable struct Liquid
    Hₗ::Float64
    ρₗ::Float64
    Cpₗ::Float64
    αₗ::Float64
    μₗ::Float64
    σ::Float64
    ζ_L::Array{Float64,1}
    ζ_R::Array{Float64,1}
    Q_B::Array{Float64,1}
    ξarrays::Array{Array{Float64,1},1}
    θarrays::Array{Array{Float64,1},1}
end

"""
Vapor is a struct containing vapor properties at a ref temperature
    Hᵥ::Float64             heat transfer coefficient between the wall and the pure vapor bubble
    k::Float64              heat conductivity of liquid in the film
    δmin::Float64           the delta with maximum heat transfer coefficient in the H interpolation
    P::Array{Float64,1}     pressure in each vapor
    Af_L, Af_R::Vector{Float64}  cross-sectional film areas (m²); use getδFromδarea for thickness
end
"""
@with_kw mutable struct Vapor
    ad_fac::Float64                 = 1.1
    Hᵥ::Float64                     = 0.0
    k::Float64
    δmin::Float64                   = 2e-6
    Eratio_plus::Float64            = 0.15
    Eratio_minus::Float64           = 0.0
    P::Array{Float64,1}
    Af_L::Array{Float64,1}
    Af_R::Array{Float64,1}
    Lfilm_start::Array{Float64,1}
    Lfilm_end::Array{Float64,1}
end

"""
Wall is a struct containing wall properties
    ΔTthres::Float64                superheat threshold to trigger boiling
    Xstations::Array{Float64,1}     locations of boiling stations on the wall
    Xarray::Array{Float64,1}        wall discrete location points for immersed boundary method
    θarray::Array{Float64,1}        wall discrete temperature points for immersed boundary method
end
"""
@with_kw mutable struct Wall
    fluid_type::String
    boil_type::String
    power::Float64
    boil_interval::Float64     = 1.0
    Rn::Float64                = 3e-6
    L_newbubble::Float64
    Xstations::Array{Float64,1}
    boiltime_stations::Array{Float64,1}
    Xarray::Array{Float64,1}
    θarray::Array{Float64,1}
    curvarray::Array{Float64,1}
end

"""
Mapping is a struct containing interpolation data
    θ_interp_walltoliquid   temperature interpolation from wall to OHP
    curv_interp_walltoliquid tube curvature interpolation from wall to OHP
    θ_interp_liquidtowall   temperature interpolation from OHP to wall
    H_interp_liquidtowall   heat transfer coefficient interpolation from OHP to wall
    P_interp_liquidtowall   pressure interpolation from OHP to wall
"""
mutable struct Mapping
    θ_interp_walltoliquid    :: AbstractInterpolation
    curv_interp_walltoliquid :: AbstractInterpolation
    θ_interp_liquidtowall    :: AbstractInterpolation
    H_interp_liquidtowall    :: AbstractInterpolation
    P_interp_liquidtowall    :: AbstractInterpolation
    heightg_interp           :: AbstractInterpolation
end

"""
Cache is a struct containing
    boil_hist    ::Vector{Any}
"""
mutable struct Cache
    boil_hist    ::Vector{Any}
    mass ::Float64
end

abstract type AbstractPHP end

"""
UniformPhpSystem is a struct containing
    tube    ::UniformTube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    mapping ::Mapping
"""
mutable struct UniformPhpSystem <: AbstractPHP
    tube    ::UniformTube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    propconvert :: PropConvert
    mapping ::Mapping
    cache   ::Cache
end

"""
This oscillating heat pipe system supports area changes along the ξ direction.
"""
mutable struct VariablePhpSystem <: AbstractPHP
    tube    ::VariableTube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    propconvert :: PropConvert
    mapping ::Mapping
    cache   ::Cache
end

"""
PhpSystem_nomapping is a struct containing
    UniformTube or VariableTube
    Liquid
    Vapor
    Wall
It is used while initializing either PHP system before its Mapping is available.
"""
mutable struct PhpSystem_nomapping{T<:Union{UniformTube,VariableTube}} <: AbstractPHP
    tube    ::T
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    propconvert :: PropConvert
end

function Base.show(io::IO, sys::UniformPhpSystem)
    N = sys.tube.N
    fluidtype = sys.propconvert.fluid_type
    println(io, "$N point OHP system filled with $fluidtype")
end
