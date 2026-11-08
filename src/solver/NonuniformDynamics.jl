export dynamicsmodel

"""
Physics calculations for an OHP channel with non-uniform geometery

Params:
    - u :: Array{Float64,1} State vector, full details on what 
    is contained within can be found in StateVectors.jl -> XMAL_to_u
    - php_sys :: VariableAreaPhpSystem

Returns:
    - du :: Array{Float64,1} first derivative of the state vector w.r.t. time
"""
function dynamicsmodel(u::Array{Float64,1}, php_sys::VariableAreaPhpSystem)
    
    du = zeros(size(u))

    
    sys :: VariableAreaPhpSystem = php_sys

    # 1st pass will only support closed tubes
    @assert sys.tube.closedornot == true;

    n_liquid_slugs = length(sys.liquid.Xp)

    Af_new = Af_at_interfaces(sys)


end