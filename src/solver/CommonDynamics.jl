export Af_at_interfaces, liquidmodel

"""
    Calculate the area of newly deposited films at each of the liquid interfaces.
    Returns `(Af_left, Af_right)`, with one area per liquid slug on each side.
    This can vary by interface due to varying geometry and direction of travel.
"""
function Af_at_interfaces(sys::Union{UniformPhpSystem, VariablePhpSystem}, liquid_ξs_left, liquid_ξs_right, liquid_char_velocity)
    @unpack Q_B = sys.liquid
    @unpack Af_L, Af_R = sys.vapor

    δ_left, δ_right = δf_new_at_interfaces(sys, liquid_ξs_left, liquid_ξs_right, liquid_char_velocity)
    Af_deposit_left = getδarea.(Ac(sys.tube, liquid_ξs_left), Dh(sys.tube, liquid_ξs_left), δ_left)
    Af_deposit_right = getδarea.(Ac(sys.tube, liquid_ξs_right), Dh(sys.tube, liquid_ξs_right), δ_right)

    Af_left = zeros(length(Q_B))
    Af_right = zeros(length(Q_B))
    indices = sys.tube.closedornot ? eachindex(Q_B) : (2:length(Q_B)-1)
    # ordering of slugs and bubbles in sys.liquid depends on whether our tube is open or closed
    left_vapor = sys.tube.closedornot ? (i -> i) : (i -> i - 1)
    right_vapor = sys.tube.closedornot ? (i -> mod1(i + 1, length(Q_B))) : (i -> i)

    for i in indices
        # Use the left vapor's "end" film if the interface is not depositing.
        Af_left[i] = Q_B[i] > 0 ? Af_deposit_left[i] : Af_R[left_vapor(i)]

        # Use the right vapor's "start" film if the interface is not depositing.
        Af_right[i] = Q_B[i] < 0 ? Af_deposit_right[i] : Af_L[right_vapor(i)]
    end

    return Af_left, Af_right
end

function δf_new_at_interfaces(sys::Union{UniformPhpSystem, VariablePhpSystem}, liquid_ξs_left, liquid_ξs_right, liquid_char_velocity)
    @unpack μₗ, σ = sys.liquid
    mean_char_velocity = mean(abs.(liquid_char_velocity))
    Ca_mean = getCa.(μₗ, σ, mean_char_velocity)
    δs_left = Catoδ.(Dh(sys.tube, liquid_ξs_left), Ca_mean; adjust_factor=sys.vapor.ad_fac)
    δs_right = Catoδ.(Dh(sys.tube, liquid_ξs_right), Ca_mean; adjust_factor=sys.vapor.ad_fac)
    return δs_left, δs_right
end

"""
    liquidmodel(sys::Union{UniformPhpSystem, VariablePhpSystem})

Model conduction & convection within liquid slugs (equation (9), Yuxuan 2024).

Physical assumptions:
- Internal temperature points are evenly spaced
- Convection coefficient is the same constant for all slugs
"""
function liquidmodel(sys::Union{UniformPhpSystem, VariablePhpSystem})
    @unpack θarrays,ξarrays,Hₗ,αₗ,Cpₗ,ρₗ = sys.liquid
    @unpack θ_interp_walltoliquid = sys.mapping
    Lt = _tube_length(sys.tube)

    du = 0 .* θarrays

    kₗ = αₗ * ρₗ * Cpₗ

    for i in eachindex(θarrays)
        ξs = ξarrays[i]
        dξs = mod.(diff(ξs), Lt)
        Δξ = sum(dξs) / length(dξs)
        
        # TODO remove once the 
        # meshing of non-uniform nodes is done
        _verify_spacing_is_uniform(dξs, Lt)

        ΔΘs = map(θ_interp_walltoliquid, ξs) - θarrays[i]

        areas = Ac(sys.tube, ξs)
        perimeters = Pc(sys.tube, ξs)

        LHS_unsteady_term_scaling = ρₗ .* Cpₗ .* areas
        RHS_conduction = kₗ .* areas .* laplacian(θarrays[i]) ./ Δξ ./ Δξ
        RHS_convection = perimeters .* Hₗ .* ΔΘs

        du[i] = (RHS_conduction + RHS_convection) ./ LHS_unsteady_term_scaling
    end

    return du
end

function _verify_spacing_is_uniform(dξs, Lt::Real)
    dξm = sum(dξs) / length(dξs)
    all(Δξ -> isfinite(Δξ) && Δξ > 0 &&
        isapprox(Δξ, dξm; rtol=sqrt(eps(Float64)), atol=16eps(Float64)*Lt), dξs) ||
        throw(ArgumentError("liquid slugs require uniformly spaced spatial nodes"))
    return nothing
end

"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""
function laplacian(u::Vector{Float64})
    unew = zeros(size(u))

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))
    d[1] = -1
    d[end]= -1

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    return (unew)
end
