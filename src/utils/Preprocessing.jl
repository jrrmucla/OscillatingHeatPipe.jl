export SimulationResult,randomζ,initialize_ohpsys,store!

"""
    randomζ(tube; numofslugs=DEFAULT_SLUGNUM,
            chargeratio=DEFAULT_LIQUID_CHARGE_RATIO, σ_charge=DEFAULT_SIGMA_CHARGE)
    randomζ(V, Vmin, closedornot; kwargs...)

Randomly generate liquid slugs in volume coordinate ζ

Params:
- tube, either a UniformTube or VariableTube
- chargeratio, expected volume fraction of liquid in the OHP
- σ_charge, standard deviation in charge ratio  

Returns:
- ζ0_L, ζ0_R, Q_B, actual_charge_ratio
"""
function randomζ(tube::UniformTube; kwargs...)
    return randomζ(tube.Ac * tube.L, tube.Ac * tube.Dh, tube.closedornot; kwargs...)
end

function randomζ(tube::VariableTube; kwargs...)
    V = map_ξ_to_ζ(tube, _tube_length(tube))
    Vmin = minimum(tube.Ac_geometric_knots .* tube.Dh_geometric_knots)
    return randomζ(V, Vmin, tube.closedornot; kwargs...)
end

function randomζ(V::Real, Vmin::Real, closedornot::Bool;
                 numofslugs=DEFAULT_SLUGNUM,
                 chargeratio=DEFAULT_LIQUID_CHARGE_RATIO,
                 σ_charge=DEFAULT_SIGMA_CHARGE)
    V > 0 || throw(ArgumentError("Total tube volume must be positive"))
    Vmin >= 0 || throw(ArgumentError("Minimum slug volume must be nonnegative"))
    numofslugs > 0 || throw(ArgumentError("Number of slugs must be positive"))
    (closedornot || numofslugs > 1) ||
        throw(ArgumentError("An open tube requires at least two boundary slugs"))

    σ_per_section = σ_charge * V / sqrt(numofslugs)
    V_per_section = V / numofslugs
    V_per_slug = V_per_section * chargeratio
    
    volumes = abs.((rand(numofslugs) .- 0.5) .* sqrt(12) .* σ_per_section .+ V_per_slug)
    minimum(volumes) > Vmin && maximum(volumes) < V_per_section ||
        error("Generation of random slugs failed.")
        # either we are trying to generate too-small slugs
        # or slugs so large they'd overlap with each other

    # place one slug at the start of each section
    ζ_L = collect(0:numofslugs-1) .* V_per_section
    ζ_R = ζ_L .+ volumes

    if closedornot
        # shift them randomly down the tube, then apply PCB
        displacement = V * rand()
        ζ_L .= mod.(ζ_L .+ displacement, V)
        ζ_R .= mod.(ζ_R .+ displacement, V)
    else
        # place the final slug against the far wall
        ζ_R[end] = V
        ζ_L[end] = V - volumes[end]
    end

    Q_B = zeros(numofslugs) # no initial velocity
    actual_charge_ratio = sum(volumes) / V
    return ζ_L, ζ_R, Q_B, actual_charge_ratio
end

"""
give d and shape of tube, return perimeter and cross-section area
"""
function peri_Ac(d::Float64,tubeshape::String)
    if tubeshape=="square"
        peri = d*4
        Ac = d*d
    elseif tubeshape=="circle"
        peri = d*π
        Ac = d*d*π/4
    end

    peri,Ac
end

"""
    initialize_ohpsys(sys::ILMSystem,p_fluid,power;[closedornot=true,
                                                    boil_waiting_time=1.0,
                                                    Rn_boil=3e-6,
                                                    inertia_f=1.3,
                                                    d=1e-3,
                                                    tube=nothing,
                                                    tubeshape="square",
                                                    Nu=3.6,
                                                    slugnum=30,
                                                    δfilm_relative=0.04,
                                                    film_fraction=0.3,
                                                    g = [0,0],
                                                    ηplus=0.6,ηminus=0,
                                                    nucleatenum=250,
                                                    L_newbubble=6e-3,
                                                    ch_ratio=0.42,
                                                    σcharge=0.01]) -> UniformPhpSystem or VariablePhpSystem

Initialize a `UniformPhpSystem`, or a `VariablePhpSystem` when given a
`VariableTube` through `tube=...`. Without a tube, construct the usual uniform
tube from `d`, `tubeshape`, `g`, and `closedornot`. A supplied tube owns these
properties; its point count and length must match the forcing shape. Variable
geometry knots must align with forcing segment endpoints. Tubes do not store
Cartesian coordinates, so this verifies arc-coordinate alignment only.
For variable tubes, the first hydraulic diameter sets the scalar liquid
heat-transfer coefficient; initial film thickness uses the local diameter.
Initializes the tube, liquid slugs, vapor regions, wall, and the mappings.
The inputs are `sys` the embedding substrate of type `ILMSystem` (which is expected to contain the
OHP geometry as a line source as the last element in the heating models array),
a tuple of fluid properties `p_fluid`, and the heater power input `power`.
"""
function initialize_ohpsys(sys::ILMSystem,p_fluid,power;closedornot=DEFAULT_CLOSED,
                                                        boil_waiting_time=DEFAULT_BOIL_WAITING_TIME,
                                                        Rn_boil=DEFAULT_RN_BOIL,
                                                        inertia_f=DEFAULT_INERTIA_F,
                                                        d=DEFAULT_D,
                                                        tube=nothing,
                                                        tubeshape=DEFAULT_TUBESHAPE,
                                                        Nu=DEFAULT_NU,
                                                        slugnum=DEFAULT_SLUGNUM,
                                                        δfilm_relative=DEFAULT_DFILM_RELATIVE,
                                                        film_fraction=DEFAULT_FILM_FRACTION,
                                                        g = DEFAULT_G,
                                                        g_potential_fn=getgh,
                                                        ηplus=DEFAULT_ETAPLUS,
                                                        ηminus=DEFAULT_ETAMINUS,
                                                        nucleatenum = DEFAULT_NUCLEATENUM,
                                                        L_newbubble = DEFAULT_L_NEWBUBBLE,
                                                        ch_ratio=DEFAULT_LIQUID_CHARGE_RATIO,
                                                        σcharge=DEFAULT_SIGMA_CHARGE)

    # Validate geometry before constructing fluid-property interpolations.
    ohp = _get_ohp_from_forcing_list(sys)
    L = arclength(ohp.shape)
    N = numpts(ohp.shape)
    if tube === nothing
        peri, Ac = peri_Ac(d, tubeshape)
        tube = UniformTube(d, peri, Ac, L, g, closedornot, N)
    else
        _verify_tube_forcing_alignment(tube, ohp.shape)
        g = tube.g
        closedornot = tube.closedornot
        d = tube isa UniformTube ? tube.Dh : first(tube.Dh_geometric_knots)
    end
    L = _tube_length(tube)

    # unpack CoolProp Properties
    @unpack fluid_type,Tref,kₗ,ρₗ,Cpₗ,αₗ,μₗ,σ = p_fluid  

    # PropConvert
    # an interpolation between different working fluid properties (currently assuming saturated gas curve)
    propconvert = PropConvert(fluid_type)

    # Liquid
    # TODO: This is fine for uniform tube but must be allowed to vary for a uniform tube
    # plus shouldn't be calculate Nu based on slug properties?
    Hₗ = p_fluid.kₗ/d * Nu # Nusselt number given

    ζ0_L,ζ0_R,Q0_B,liquid_realratio = randomζ(tube,numofslugs=slugnum,chargeratio=ch_ratio,σ_charge=σcharge) # chargeratio here is nominal volume fraction for liquid slugs only, it is slightly different from actual volume fraction.
    ξ0_L = map_ζ_to_ξ(tube, ζ0_L)
    ξ0_R = map_ζ_to_ξ(tube, ζ0_R)
    ξ0 = collect(zip(ξ0_L, ξ0_R))
    Xarrays,θarrays = constructXarrays(ξ0,N,Tref,L)
    
    liquids=Liquid(Hₗ,ρₗ,Cpₗ,αₗ,μₗ,σ,ζ0_L,ζ0_R,Q0_B,Xarrays,θarrays)

    # Vapor
    @unpack TtoP = propconvert

    Lvaporplug = ξ_to_vapor_length(ξ0_L, ξ0_R,L,tube.closedornot)
    P_initial = zero(Lvaporplug) .+ TtoP(Tref)
    Xpvapor = getXpvapor(ξ0, tube.closedornot)
    Af_L_initial, Af_R_initial = _initial_film_areas(tube, Xpvapor, δfilm_relative)
    Lfilm_start_initial = 0.5 .* film_fraction .* Lvaporplug
    Lfilm_end_initial   = deepcopy(Lfilm_start_initial)
    vapors=Vapor(ad_fac=inertia_f,k = p_fluid.kₗ,P=P_initial,Af_L=Af_L_initial,Af_R=Af_R_initial,Lfilm_start=Lfilm_start_initial,Lfilm_end=Lfilm_end_initial,Eratio_plus=ηplus,Eratio_minus=ηminus);

    # Wall

    # boiling stations
    ξ_boiling_stations = sort(rand(nucleatenum) .* L);
    ξ_station_time = zeros(nucleatenum);
    boil_type = "wall T"
    boil_interval = boil_waiting_time
    
    # thermal points of interaction with plate system
    Xwallarray,θwallarray,curvwallarray = constructwallXθarray(arccoordmid(ohp.shape),Tref,curvature(ohp.shape));
    wall = Wall(boil_interval=boil_interval,fluid_type=fluid_type,boil_type=boil_type,power=power,L_newbubble=L_newbubble,Xstations=ξ_boiling_stations,boiltime_stations=ξ_station_time,Xarray=Xwallarray,θarray=θwallarray,curvarray=curvwallarray,Rn=Rn_boil);

    # Mapping
    @unpack x,y = ohp.transform(ohp.shape)
    sys0_nomapping = PhpSystem_nomapping(tube,liquids,vapors,wall,propconvert);
    
    # Gravitational Potential (as a function of space and time)
    # we require the mapping function's signiture ~ fn(g, x, y)
    @assert applicable(g_potential_fn, g, x, y)
    ht = g_potential_fn(g,x,y)
    heightg_interp = linear_interpolation(Xwallarray,ht,extrapolation_bc = Line())
    
    mapping = initialize_interpolations(sys0_nomapping, heightg_interp)

    # Cache
    Mfilm_left, Mfilm_right = getMfilm(sys0_nomapping)
    totalmass = sum(getMvapor(sys0_nomapping)) + sum(getMliquid(sys0_nomapping)) + sum(Mfilm_left .+ Mfilm_right)
    boil_hist_int = []
    cache = Cache(boil_hist_int,totalmass)

    # UniformPhpSystem or VariablePhpSystem
    system_constructor = tube isa VariableTube ? VariablePhpSystem : UniformPhpSystem
    sys0 = system_constructor(tube,liquids,vapors,wall,propconvert,mapping,cache);

    sys0
end

function _verify_tube_forcing_alignment(tube, shape)
    tube isa Union{UniformTube,VariableTube} ||
        throw(ArgumentError("tube must be a UniformTube or VariableTube"))
    tube.N == numpts(shape) ||
        throw(ArgumentError("tube point count does not match the OHP forcing shape"))
    if tube isa VariableTube
        length(tube.ξ_geometric_knots) >= 2 ||
            throw(ArgumentError("tube geometry requires at least two knots"))
    end
    isapprox(_tube_length(tube), arclength(shape); rtol=1e-8, atol=1e-12) ||
        throw(ArgumentError("tube length does not match the OHP forcing shape"))
    tube isa UniformTube && return nothing

    knots = tube.ξ_geometric_knots
    length(knots) >= 2 && first(knots) == 0 && all(diff(knots) .> 0) ||
        throw(ArgumentError("tube geometry knots must increase strictly from zero"))
    for values in (tube.Dh_geometric_knots, tube.Pc_geometric_knots, tube.Ac_geometric_knots)
        length(values) == length(knots) && all(value -> isfinite(value) && value > 0, values) ||
            throw(ArgumentError("each geometry knot must have positive finite Dh, Pc, and Ac values"))
    end
    forcing_endpoints = [0.0; arccoord(shape)]
    tolerance = max(1e-12, 1e-8 * arclength(shape))
    for knot in knots
        index = searchsortedfirst(forcing_endpoints, knot)
        candidates = max(1, index - 1):min(length(forcing_endpoints), index)
        any(abs(forcing_endpoints[i] - knot) <= tolerance for i in candidates) ||
            throw(ArgumentError("tube geometry knot at $knot does not align with a forcing segment endpoint"))
    end
    return nothing
end

function _initial_film_areas(tube::UniformTube, vapor_interfaces, relative_thickness)
    film_area = getδarea(tube.Ac, tube.Dh, relative_thickness * tube.Dh / 2)
    return fill(film_area, length(vapor_interfaces)), fill(film_area, length(vapor_interfaces))
end

function _initial_film_areas(tube::VariableTube, vapor_interfaces, relative_thickness)
    diameter, _, area = construct_geometric_interpolations(tube)
    film_area(position) = getδarea(area(position), diameter(position), relative_thickness * diameter(position) / 2)
    return film_area.(first.(vapor_interfaces)), film_area.(last.(vapor_interfaces))
end

function _get_ohp_from_forcing_list(sys::ILMSystem)
    ohp = sys.forcing["heating models"][end] # should be at end of list
    ohp isa LineForcingModel || error("OHP not in expected place in forcing list")
    return ohp
end


mutable struct SimulationResult
    tube_hist_t      ::Vector{Any}
    tube_hist_u      ::Vector{Any}
    tube_hist_θwall  ::Vector{Any}
    boil_hist        ::Vector{Any}
    plate_T_hist     ::Vector{Any}
    integrator_tube  ::Any
    integrator_plate ::Any
    integrator_tube_resume  ::Any
    integrator_plate_resume ::Any    
    # grid             ::Any
end

function SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate)
    integrator_tube_resume = deepcopy(integrator_tube)
    integrator_plate_resume = deepcopy(integrator_plate)
    SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,integrator_tube_resume,integrator_plate_resume)
    # SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,integrator_tube_resume,integrator_plate_resume,grid)
end

function SimulationResult(int_tube,int_plate)
    boil_hist = []
    plate_T_hist = []
    tube_hist_u  = []
    tube_hist_t = []
    tube_hist_θwall = []
    integrator_tube = deepcopy(int_tube)
    integrator_plate = deepcopy(int_plate)

    return SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate)
end

function store!(sr,integrator_tube,integrator_plate)
        
        append!(sr.boil_hist,deepcopy(integrator_tube.p.cache.boil_hist));
        integrator_tube.p.cache.boil_hist = []

        push!(sr.plate_T_hist,deepcopy(temperature(integrator_plate)));
        push!(sr.tube_hist_θwall,deepcopy(integrator_tube.p.wall.θarray))
        push!(sr.tube_hist_u,deepcopy(integrator_tube.u));
        push!(sr.tube_hist_t,deepcopy(integrator_tube.t));
        sr.integrator_tube_resume = deepcopy(integrator_tube)
        sr.integrator_plate_resume = deepcopy(integrator_plate)

end
