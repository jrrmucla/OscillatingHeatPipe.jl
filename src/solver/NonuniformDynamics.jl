export dynamicsmodel

"""
Physics calculations for an OHP channel with non-uniform geometery

Params:
    - u :: Array{Float64,1} state vector
            [ζ_L (n); ζ_R (n); Q_B (n); M (n*); Af_L (n*); Af_R (n*); Lfilm_start (n*); Lfilm_end (n*)]
             m^3      m^3      m^3/s    kg      m^2        m^2        m                 m

    - php_sys :: VariablePhpSystem

Returns:
    - du :: Array{Float64,1} first derivative of the state vector w.r.t. time
            [Q_L (n); Q_R (n); dQᵦdt (n); ṁ (n*); Ȧf_L;  Ȧf_R;  L̇f_L; L̇f_R]
             m^3/s    m^3/s    m^3/s^2    kg/s    m^2/s  m^2/s  m/s   m/s
"""
function dynamicsmodel(u::Array{Float64,1}, php_sys::VariablePhpSystem)
    sys::VariablePhpSystem = php_sys
    # 1st pass will only support closed tubes
    @assert sys.tube.closedornot == true;

    du = zeros(size(u))
    
    @unpack closedornot = sys.tube
    @unpack σ,μₗ,ρₗ,ζ_L,ζ_R,Q_B = sys.liquid
    @unpack P,Eratio_plus,Eratio_minus,Af_L,Af_R,Lfilm_start,Lfilm_end,ad_fac = sys.vapor
    @unpack L_newbubble = sys.wall

    Vt::Float64 = _tube_volume(sys.tube)
    Lt::Float64 = _tube_length(sys.tube)

    # cache some information about our liquid slugs
    n_liquid_slugs::Int = length(sys.liquid.ζ_L)
    liquid_ξ_left = map_ζ_to_ξ.(Ref(sys.tube), ζ_L)
    liquid_ξ_right = map_ζ_to_ξ.(Ref(sys.tube), ζ_R)
    liquid_lengths = ξ_to_liquid_length(liquid_ξ_left, liquid_ξ_right, Lt)
    liquid_volumes = ζ_to_liquid_volume(ζ_L, ζ_R, Vt)

    # cache some information about the vapor bubbles
    n_vapor_bubbles::Int = length(sys.vapor.Af_L)
    vapor_lengths = ξ_to_vapor_length(liquid_ξ_left, liquid_ξ_right, Lt, closedornot)

    # Film Deposition
    Ac_left_interface::Vector{Float64} = Ac(sys.tube, liquid_ξ_left)
    Ac_right_interface::Vector{Float64} = Ac(sys.tube, liquid_ξ_right)
    old_velocity_left_interface = Q_B ./ Ac_left_interface
    old_velocity_right_interface = Q_B ./ Ac_right_interface
    characteristic_velocity_film_deposition = (old_velocity_left_interface .+ old_velocity_right_interface) ./ 2.0
    ## this value gets reduced to a common Ca for all slugs

    ## Yuxuan 2024, equation (3)
    Af_left_interface::Vector{Float64}, Af_right_interface::Vector{Float64} = Af_at_interfaces(
        sys, liquid_ξ_left, liquid_ξ_right, characteristic_velocity_film_deposition)

    # Conservation of Mass
    Q_left_interface::Vector{Float64} = Q_B ./ (1.0 .- ((Af_left_interface ./ Ac_left_interface)))
    Q_right_interface::Vector{Float64} = Q_B ./ (1.0 .- ((Af_right_interface ./ Ac_right_interface)))

    # Q = dζ/dt
    du[               1 :   n_liquid_slugs] = Q_left_interface
    du[n_liquid_slugs+1 : 2*n_liquid_slugs] = Q_right_interface

    # Conservation of Momentum
    ##
    ##
    ##

    return du
end
