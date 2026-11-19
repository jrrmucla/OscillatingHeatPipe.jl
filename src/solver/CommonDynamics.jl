export Af_at_interfaces

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
