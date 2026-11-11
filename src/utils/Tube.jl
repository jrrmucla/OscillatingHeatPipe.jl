export UniformTube

"""
UniformTube is a struct containing tube geometries
    Dh          ::Float64   hydraulic diameter
    Pc          ::Float64   channel perimeter (may be an arbitrary shape so cannot be derived from d)
    Ac          ::Float64   cross-sectional area
    L           ::Float64   tube one dimensional length
    g           ::Float64   gravity vector along the plate
    closedornot ::Bool      if the tube is closed loop or not (open loop)
"""
struct UniformTube
    Dh::Float64
    Pc::Float64
    Ac::Float64
    L::Float64
    g::Vector{Float64}
    closedornot::Bool
    N::Int64
end

"""
Nonuniform tube and its immutable cached lookups. Constructor function
copies geometry inputs to builds all interpolants and coordinate maps. 
Treat the stored knots as read-only; construct a new tube to change its geometry.
"""
struct VariableTube{A,B,C,D,E}
    ξ_geometric_knots::Array{Float64,1}
    Dh_geometric_knots::Array{Float64,1}
    Pc_geometric_knots::Array{Float64,1}
    Ac_geometric_knots::Array{Float64,1}
    g::Vector{Float64}
    closedornot::Bool
    N::Int64
    map_ξ_to_Dh::A
    map_ξ_to_Pc::B
    map_ξ_to_Ac::C
    map_ξ_to_ζ::D
    map_ζ_to_ξ::E
end

# Constructor for VariableTube
# we don't want to manually set the internal geometric maps,
# so they are computed automatically given the parameters in this function
# soon we won't even want users to set all fields 1-by-1, instead give a shape
# circular/rectangular, and D / W array, H array
function VariableTube(ξknots, diameters, perimeters, areas, g, closedornot, N)
    ξknots = collect(Float64, ξknots)
    diameters = collect(Float64, diameters)
    perimeters = collect(Float64, perimeters)
    areas = collect(Float64, areas)
    length(diameters) == length(perimeters) == length(areas) == length(ξknots) ||
        throw(ArgumentError("all geometry knot arrays must have equal lengths"))
    forward, inverse = construct_exact_coordinate_maps(ξknots, areas)
    diameter_map = linear_interpolation(ξknots, diameters)
    perimeter_map = linear_interpolation(ξknots, perimeters)
    area_map = linear_interpolation(ξknots, areas)
    return VariableTube(ξknots, diameters, perimeters, areas,
                        collect(Float64, g), Bool(closedornot), Int64(N),
                        diameter_map, perimeter_map, area_map, forward, inverse)
end

_tube_length(tube::UniformTube) = tube.L
_tube_length(tube::VariableTube) = last(tube.ξ_geometric_knots)

Dh(tube::UniformTube, ξ::Real) = tube.Dh
Pc(tube::UniformTube, ξ::Real) = tube.Pc
Ac(tube::UniformTube, ξ::Real) = tube.Ac
Dh(tube::VariableTube, ξ::Real) = tube.map_ξ_to_Dh(ξ)
Pc(tube::VariableTube, ξ::Real) = tube.map_ξ_to_Pc(ξ)
Ac(tube::VariableTube, ξ::Real) = tube.map_ξ_to_Ac(ξ)

map_ξ_to_ζ(tube::UniformTube, ξ::Real) = tube.Ac * ξ
map_ζ_to_ξ(tube::UniformTube, ζ::Real) = ζ / tube.Ac
map_ξ_to_ζ(tube::VariableTube, ξ::Real) = tube.map_ξ_to_ζ(ξ)
map_ζ_to_ξ(tube::VariableTube, ζ::Real) = tube.map_ζ_to_ξ(ζ)

for lookup in (:Dh, :Pc, :Ac, :map_ξ_to_ζ, :map_ζ_to_ξ)
    @eval $lookup(tube::Union{UniformTube,VariableTube}, positions::AbstractArray) =
        $lookup.(Ref(tube), positions)
end

# Retain existing geometry consumers while returning the tube's cached lookups.
construct_geometric_interpolations(tube::VariableTube) =
    (tube.map_ξ_to_Dh, tube.map_ξ_to_Pc, tube.map_ξ_to_Ac)
construct_geometric_interpolations(tube::UniformTube) =
    (ξ -> Dh(tube, ξ), ξ -> Pc(tube, ξ), ξ -> Ac(tube, ξ))
construct_exact_coordinate_maps(tube::VariableTube) =
    (tube.map_ξ_to_ζ, tube.map_ζ_to_ξ)
construct_exact_coordinate_maps(tube::UniformTube) =
    (ξ -> map_ξ_to_ζ(tube, ξ), ζ -> map_ζ_to_ξ(tube, ζ))

_tube_diameters(tube::Union{UniformTube,VariableTube}, positions) = Dh(tube, positions)
_tube_perimeters(tube::Union{UniformTube,VariableTube}, positions) = Pc(tube, positions)
_tube_areas(tube::Union{UniformTube,VariableTube}, positions) = Ac(tube, positions)

"""
    construct_exact_coordinate_maps(ξknots::AbstractVector, areas::AbstractVector)

Return `(map_ξ_to_ζ, map_ζ_to_ξ)` for `ζ = ∫_0^ξ A(s) ds`.
Used by the tube constructor; passing a tube instead returns its cached maps.
Both maps accept scalars or arrays within the tube. Knots must start at zero
and increase strictly, with positive areas.

One-point Gauss quadrature is exact for each linear
area segment, including partial segments, up to floating-point roundoff.

Use broadcasting into a preallocated array (`dest .= map_ξ_to_ζ.(src)`) to reuse output storage.
"""
function construct_exact_coordinate_maps(ξknots::AbstractVector, areas::AbstractVector)
    ξknots = collect(Float64, ξknots)
    areas = collect(Float64, areas)

    # io validation
    length(ξknots) == length(areas) >= 2 ||
        throw(ArgumentError("coordinate and area knots must have equal lengths of at least two"))
    all(isfinite, ξknots) && first(ξknots) == 0 && all(>(0), diff(ξknots)) ||
        throw(ArgumentError("coordinate knots must be finite, start at zero, and increase strictly"))
    all(a -> isfinite(a) && a > 0, areas) ||
        throw(ArgumentError("area knots must be finite and positive"))

    widths = diff(ξknots)
    slopes = diff(areas) ./ widths

    # Gauss node 0 and weight 2 on [-1, 1], mapped to each segment.
    ζknots = [0.0; cumsum(widths .* (areas[1:end-1] ./ 2 .+ areas[2:end] ./ 2))]
    halfslopes = slopes ./ 2
    squared_areas = areas .^ 2
    inverse_volumes = 1 ./ diff(ζknots)
    nsegments = length(slopes)
    ξmax = last(ξknots)
    ζmax = last(ζknots)

    @inline function map_ξ_to_ζ(ξ::Real)
        0 <= ξ <= ξmax ||
            throw(DomainError(ξ, "length coordinate is outside the tube"))
        i = min(searchsortedlast(ξknots, ξ), nsegments)
        @inbounds begin
            s = ξ - ξknots[i]
            return ζknots[i] + s * (areas[i] + halfslopes[i] * s)
        end
    end

    @inline function map_ζ_to_ξ(ζ::Real)
        0 <= ζ <= ζmax ||
            throw(DomainError(ζ, "volume coordinate is outside the tube"))
        i = min(searchsortedlast(ζknots, ζ), nsegments)
        @inbounds begin
            Δζ = ζ - ζknots[i]
            # This equals Aᵢ² + 2mᵢΔζ, avoiding subtraction for decreasing area.
            fraction = Δζ * inverse_volumes[i]
            area_squared = (1 - fraction) * squared_areas[i] + fraction * squared_areas[i+1]
            s = 2Δζ / (areas[i] + sqrt(area_squared))
            return ξknots[i] + s
        end
    end

    map_ξ_to_ζ(ξ::AbstractArray) = map_ξ_to_ζ.(ξ)
    map_ζ_to_ξ(ζ::AbstractArray) = map_ζ_to_ξ.(ζ)
    return map_ξ_to_ζ, map_ζ_to_ξ
end
