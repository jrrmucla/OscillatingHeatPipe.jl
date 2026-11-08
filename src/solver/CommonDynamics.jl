export Af_at_interfaces

"""
    Calculate the area of newly deposited films at each of the liquid interfaces.
    This can vary by interface due to varying geometry and direction of travel.
"""
function Af_at_interfaces(sys::Union{UniformPhpSystem, VariableAreaPhpSystem})
    @unpack Xp, dXdt = sys.liquid
    @unpack Afilm_start, Afilm_end = sys.vapor

    δ_deposit = δf_new_at_interfaces(sys)
    ξl = collect(Iterators.flatten(Xp))
    Af_deposit = δ_deposit .* Pc(sys, ξl) 
    # this is a flattened array of Af for deposition

    next_index = circshift(1:length(Xp),-1)
    actual_Af_int = fill((0.0, 0.0), length(dXdt))

    indices = sys.tube.closedornot ? eachindex(actual_Af_int) : (2:length(actual_Af_int)-1)

    # ordering of slugs and bubbles in sys.liquid depends on whether our tube is open or closed
    left_vapor = sys.tube.closedornot ? (i -> i) : (i -> i - 1)
    right_vapor = sys.tube.closedornot ? (i -> next_index[i]) : (i -> i)

    for i in indices
        # Use the left vapor's "end" film if the interface is not depositing.
        Adeposit_left = dXdt[i][1] > 0 ? Af_deposit[2i-1] : Afilm_end[left_vapor(i)]

        # Use the right vapor's "start" film if the interface is not depositing.
        Adeposit_right = dXdt[i][end] < 0 ? Af_deposit[2i] : Afilm_start[right_vapor(i)]

        actual_Af_int[i]  =   (Adeposit_left, Adeposit_right)
    end

    return actual_Af_int
end

function δf_new_at_interfaces(sys::Union{UniformPhpSystem, VariableAreaPhpSystem})
    @unpack Xp, dXdt, μₗ, σ = sys.liquid
    V = [mean(elem) for elem in dXdt]
    V_mean = mean(abs.(V))
    Ca_mean = getCa.(μₗ, σ, V_mean)
    ξl = collect(Iterators.flatten(Xp))
    diameters = Dh(sys, ξl)
    δs = Catoδ(diameters, Ca_mean, adjust_factor=sys.vapor.ad_fac)
    return δs
end
