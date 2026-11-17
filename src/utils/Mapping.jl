export initialize_interpolations, refresh_interpolations, slug_interp

"""
    initialize_interpolations(sys::AbstractPHP, map_ξ_to_heightg::AbstractInterpolation)

Returns a Mapping for sys. Creates wall and fluid interpolations and stores
the gravitational potential interpolant. Geometry lookups live in the tube.

After initialization, use `refresh_interpolations(sys)` to reuse the
fixed tube geometry instead of constructing it again.
"""
function initialize_interpolations(sys::AbstractPHP, map_ξ_to_heightg::AbstractInterpolation)
    return assemble_interpolations(sys, map_ξ_to_heightg)
end

"""
    refresh_interpolations(sys::AbstractPHP)

Return a new Mapping using the current wall and fluid state. Rebuilds wall and fluid interps. 
Does not rebuild geometric interps or the gravity-potential interp.
"""
function refresh_interpolations(sys::AbstractPHP)
    mapping = sys.mapping
    return assemble_interpolations(sys, mapping.heightg_interp)
end

function assemble_interpolations(sys::AbstractPHP, map_ξ_to_heightg::AbstractInterpolation)
    map_ξ_to_θw, map_ξ_to_wall_curvature, map_ξ_to_θf, map_ξ_to_Hf, map_ξ_to_Pf =
        sys_interpolation(sys)
    return Mapping(map_ξ_to_θw, map_ξ_to_wall_curvature,
                   map_ξ_to_θf, map_ξ_to_Hf, map_ξ_to_Pf,
                   map_ξ_to_heightg)
end

function sys_interpolation(sys)
    if sys.tube.closedornot == true
            return sys_interpolation_closedloop(sys)
        else
            return sys_interpolation_openloop(sys)
    end
end

function sys_interpolation_openloop(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)

    # these maps are applied once in dynamics, once here
    # likely worth caching since the computation isn't trivial
    ξ_L = map_ζ_to_ξ(sys.tube, sys.liquid.ζ_L)
    ξ_R = map_ζ_to_ξ(sys.tube, sys.liquid.ζ_R)
    ξp = collect(zip(ξ_L, ξ_R))

    @unpack PtoT = sys.propconvert
    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    Af_L = sys.vapor.Af_L
    Af_R = sys.vapor.Af_R
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end
    ξp_vapor = getXpvapor(ξp,sys.tube.closedornot)

    H_film_start = _film_heat_transfer(sys, Af_L, first.(ξp_vapor))
    H_film_end = _film_heat_transfer(sys, Af_R, last.(ξp_vapor))
    H_vapor = sys.vapor.Hᵥ
    H_liquid = sys.liquid.Hₗ

    Nvapor = length(P)

    append!(X_inner,sys.liquid.ξarrays[1])
    append!(θ_inner,sys.liquid.θarrays[1])
    append!(H_inner,H_liquid .* ones(length(sys.liquid.ξarrays[1])))
    append!(X_inner_pres,[ξp[1][1],ξp[1][end]])
    append!(P_inner,[P[1], P[1]])

    for i = 1:Nvapor

        # vapors
                append!(X_inner,[ξp_vapor[i][1],ξp_vapor[i][1]+Lfilm_start[i],
                ξp_vapor[i][1]+Lfilm_start[i],ξp_vapor[i][end]-Lfilm_end[i],
                ξp_vapor[i][end]-Lfilm_end[i],ξp_vapor[i][end]])

                append!(θ_inner,[θ[i],θ[i],θ[i],θ[i],θ[i],θ[i]])
                append!(H_inner,[H_film_start[i], H_film_start[i],
                H_vapor,H_vapor,
                H_film_end[i],H_film_end[i]])

                append!(X_inner_pres,[ξp_vapor[i][1],ξp_vapor[i][end]])
                append!(P_inner,[P[i], P[i]])
    
        # liquids
            append!(X_inner,sys.liquid.ξarrays[i+1])
            append!(θ_inner,sys.liquid.θarrays[i+1])
            append!(H_inner,H_liquid .* ones(length(sys.liquid.ξarrays[i+1])))
            append!(X_inner_pres,[ξp[i+1][1],ξp[i+1][end]])
            append!(P_inner,[P[i], P[i]])
        end

    # ``` extend wall Xarray by adding its 0.0 point (plan to write this as a seperate function)```
    extend_wall_Xarray = deepcopy(sys.wall.Xarray)
    extend_wall_θarray = deepcopy(sys.wall.θarray)

    extend_wall_Xarray = [0.0;sys.wall.Xarray;_tube_length(sys.tube)]
    extend_wall_θarray = [(sys.wall.θarray[1]+sys.wall.θarray[end])/2;sys.wall.θarray;(sys.wall.θarray[1]+sys.wall.θarray[end])/2]
    extend_wall_curvarray = [(sys.wall.curvarray[1]+sys.wall.curvarray[end])/2;sys.wall.curvarray;(sys.wall.curvarray[1]+sys.wall.curvarray[end])/2]

    Interpolations.deduplicate_knots!(X_inner,move_knots = true)
    Interpolations.deduplicate_knots!(extend_wall_Xarray,move_knots = true)
    Interpolations.deduplicate_knots!(X_inner_pres,move_knots = true)

    θ_interp_liquidtowall = linear_interpolation(X_inner, θ_inner);

    H_interp_liquidtowall = linear_interpolation(X_inner, H_inner);

    θ_interp_walltoliquid = linear_interpolation(extend_wall_Xarray, extend_wall_θarray);

    curv_interp_walltoliquid = linear_interpolation(extend_wall_Xarray, extend_wall_curvarray)

    P_interp_liquidtowall = linear_interpolation(X_inner_pres, P_inner);


    return θ_interp_walltoliquid, curv_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall

end


function sys_interpolation_closedloop(sys)
    # this function (and its openloop twin) 
    # is my op for a variety of reasons
    # count your days, sys_interpolation_closedloop.... a refactor is coming...
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)

    ξ_L = map_ζ_to_ξ(sys.tube, sys.liquid.ζ_L)
    ξ_R = map_ζ_to_ξ(sys.tube, sys.liquid.ζ_R)
    ξp = collect(zip(ξ_L, ξ_R))

    @unpack PtoT = sys.propconvert
    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    Af_L = sys.vapor.Af_L
    Af_R = sys.vapor.Af_R
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end
    ξp_vapor = getXpvapor(ξp,sys.tube.closedornot)

    H_film_start = _film_heat_transfer(sys, Af_L, first.(ξp_vapor))
    H_film_end = _film_heat_transfer(sys, Af_R, last.(ξp_vapor))
    H_vapor = sys.vapor.Hᵥ
    H_liquid = sys.liquid.Hₗ

    Nvapor = length(P)

    loop_plus_index = [2:Nvapor;1]

    max_i = argmax(ξp_vapor) # vapor plug with largest arc coordinate
    max_j = argmax(ξp_vapor[max_i]) # end of this plug with the largest coordinate (1 or 2)
  
    # Loop through the vapor and slugs (which there are equal number of)
    for i = 1:Nvapor

            if i == max_i && max_j != 2 # if this is the furthest plug and it crosses the branch cut

                X_inner_loop_append,H_inner_loop_append = XHloop_append(i,ξp_vapor[i],H_film_start[i],H_film_end[i],sys)

                append!(X_inner,X_inner_loop_append)
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,H_inner_loop_append)

                append!(X_inner_pres,[ξp_vapor[i][1], _tube_length(sys.tube), 0.0, ξp_vapor[i][end]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
            else 
                append!(X_inner,[ξp_vapor[i][1],ξp_vapor[i][1]+Lfilm_start[i],
                ξp_vapor[i][1]+Lfilm_start[i],ξp_vapor[i][end]-Lfilm_end[i],
                ξp_vapor[i][end]-Lfilm_end[i],ξp_vapor[i][end]])

                append!(θ_inner,[θ[i],θ[i],θ[i],θ[i],θ[i],θ[i]])
                append!(H_inner,[H_film_start[i], H_film_start[i],
                H_vapor,H_vapor,
                H_film_end[i],H_film_end[i]])

                append!(X_inner_pres,[ξp_vapor[i][1],ξp_vapor[i][end]])
                append!(P_inner,[P[i], P[i]])
    
        end
        

    if (sys.liquid.ξarrays[i][1] > sys.liquid.ξarrays[i][end])
        period_index = findmin(sys.liquid.ξarrays[i])[2]
        Xarrays_temp = deepcopy(sys.liquid.ξarrays[i])
        θarrays_temp = deepcopy(sys.liquid.θarrays[i])
        
        
        insert!(Xarrays_temp,period_index, _tube_length(sys.tube))
        insert!(Xarrays_temp,period_index+1, 0.0)
        insert!(θarrays_temp,period_index, sys.liquid.θarrays[i][period_index])
        insert!(θarrays_temp,period_index+1, sys.liquid.θarrays[i][period_index])
        append!(X_inner,Xarrays_temp)
        append!(θ_inner,θarrays_temp)
        
        append!(H_inner,H_liquid .* ones(length(sys.liquid.ξarrays[i]) + 2))
        

        append!(X_inner_pres,[sys.liquid.ξarrays[i][1], _tube_length(sys.tube), 0.0, sys.liquid.ξarrays[i][end]])
        P_inner_end = (_tube_length(sys.tube)-sys.liquid.ξarrays[i][1])/mod(sys.liquid.ξarrays[i][end]-sys.liquid.ξarrays[i][1],_tube_length(sys.tube)) * (P[loop_plus_index[i]]-P[i]) + P[i]
        append!(P_inner,[P[i], P_inner_end, P_inner_end, P[loop_plus_index[i]]])

        else
            append!(X_inner,sys.liquid.ξarrays[i])
            append!(θ_inner,sys.liquid.θarrays[i])
            append!(H_inner,H_liquid .* ones(length(sys.liquid.ξarrays[i])))
        end
end

        imin = argmin(X_inner)
        X_inner_final = X_inner
        θ_inner_final = θ_inner
        H_inner_final = H_inner
        if imin != 1
            X_inner_final = [X_inner[imin:end];X_inner[1:imin-1]]
            θ_inner_final = [θ_inner[imin:end];θ_inner[1:imin-1]]
            H_inner_final = [H_inner[imin:end];H_inner[1:imin-1]]
        end

        imin = argmin(X_inner_pres)
        X_inner_pres_final = X_inner_pres
        P_inner_final = P_inner
        if imin != 1
            X_inner_pres_final = [X_inner_pres[imin:end];X_inner_pres[1:imin-1]]
            P_inner_final = [P_inner[imin:end];P_inner[1:imin-1]]
        end
    
    ``` extend wall Xarray by adding its 0.0 point```
    extend_wall_Xarray = deepcopy(sys.wall.Xarray)
    extend_wall_θarray = deepcopy(sys.wall.θarray)

    extend_wall_Xarray = [0.0;sys.wall.Xarray;_tube_length(sys.tube)]
    extend_wall_θarray = [(sys.wall.θarray[1]+sys.wall.θarray[end])/2;sys.wall.θarray;(sys.wall.θarray[1]+sys.wall.θarray[end])/2]
    extend_wall_curvarray = [(sys.wall.curvarray[1]+sys.wall.curvarray[end])/2;sys.wall.curvarray;(sys.wall.curvarray[1]+sys.wall.curvarray[end])/2]


    Interpolations.deduplicate_knots!(X_inner_final,move_knots = true)
    Interpolations.deduplicate_knots!(extend_wall_Xarray,move_knots = true)
    Interpolations.deduplicate_knots!(X_inner_pres_final,move_knots = true)

    θ_interp_liquidtowall = linear_interpolation(X_inner_final, θ_inner_final);

    H_interp_liquidtowall = linear_interpolation(X_inner_final, H_inner_final);

    θ_interp_walltoliquid = linear_interpolation(extend_wall_Xarray, extend_wall_θarray);

    curv_interp_walltoliquid = linear_interpolation(extend_wall_Xarray, extend_wall_curvarray)

    P_interp_liquidtowall = linear_interpolation(X_inner_pres_final, P_inner_final);


    return θ_interp_walltoliquid, curv_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall

end

function XHloop_append(i,ξp_vapor,H_film_start,H_film_end,sys)


    Lfilm_start = sys.vapor.Lfilm_start[i]
    Lfilm_end = sys.vapor.Lfilm_end[i]

    H_vapor = sys.vapor.Hᵥ

    L = _tube_length(sys.tube)

    if ξp_vapor[2] > ξp_vapor[1]
        # println("Xp error 1")
        return error("Xp error 1")
        
    elseif mod(ξp_vapor[1],L) + Lfilm_start > L
        
        X_inner_loop_append = mod.([ξp_vapor[1],L,
                0.0,ξp_vapor[1]+Lfilm_start,
                ξp_vapor[1]+Lfilm_start,ξp_vapor[end]-Lfilm_end,
                ξp_vapor[end]-Lfilm_end,ξp_vapor[end]],L)
        X_inner_loop_append[2] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_film_start, H_film_start,
                H_vapor,H_vapor,
                H_film_end,H_film_end]
        
    elseif ξp_vapor[2] - Lfilm_end > 0.0
        
        X_inner_loop_append = mod.([ξp_vapor[1],ξp_vapor[1]+Lfilm_start,
                ξp_vapor[1]+Lfilm_start, L,
                0.0, ξp_vapor[end]-Lfilm_end,
                ξp_vapor[end]-Lfilm_end,ξp_vapor[end]],L)
        X_inner_loop_append[4] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_vapor, H_vapor,
                H_film_end,H_film_end]
        
    elseif ξp_vapor[2] > 0
        
        X_inner_loop_append = mod.([ξp_vapor[1],ξp_vapor[1]+Lfilm_start,
                ξp_vapor[1]+Lfilm_start,  ξp_vapor[end]-Lfilm_end,
                ξp_vapor[end]-Lfilm_end,L,
                0.0,ξp_vapor[end]],L)
        X_inner_loop_append[6] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_film_end,H_film_end,
                H_film_end,H_film_end]
        
    else 
        println("Xp error 2")
        return "Xp error 2"
    end

return X_inner_loop_append,H_inner_loop_append
end

function slug_interp(sys::UniformPhpSystem)
        filmLend = sys.vapor.Lfilm_end
    filmLstart = sys.vapor.Lfilm_start
    # Xpstart = [elem[1] for elem in sys.liquid.Xp]
    # Xpend = [elem[2] for elem in sys.liquid.Xp]
    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.closedornot)
    Xvaporend = [elem[2] for elem in Xpvapor]
    Xvaporstart = [elem[1] for elem in Xpvapor]
    Xfilmstart = Xvaporstart .+ filmLstart
    Xfilmend = Xvaporend .- filmLend
    Xslugs = zeros(8*size(Xvaporstart)[1])
    slug_flags = zeros(8*size(Xvaporstart)[1])
    Xfinalslugs = zeros(8*size(Xvaporstart)[1])
    slug_flags = zeros(8*size(Xvaporstart)[1])

    for i in eachindex(Xvaporstart)
        Xslugs[8i-7:8i-6] = [Xvaporstart[i];Xvaporstart[i]]
        Xslugs[8i-5:8i-4] = [Xfilmstart[i];Xfilmstart[i]]
        Xslugs[8i-3:8i-2] = [Xfilmend[i];Xfilmend[i]]
        Xslugs[8i-1:8i] = [Xvaporend[i];Xvaporend[i]]

        slug_flags[8i-7:8i-6] = [1.0; 2.0]
        slug_flags[8i-5:8i-4] = [2.0; 0.0]
        slug_flags[8i-3:8i-2] = [0.0; 2.0]
        slug_flags[8i-1:8i] =   [2.0; 1.0]

        Xslugs = mod.(Xslugs,_tube_length(sys.tube))
    end

    Xmax, imax = findmax(Xslugs)
    #     Xmin, imin = findmax(Xslugs)

    edge_status = slug_flags[imax+1]

    xtemp = splice!(Xslugs,1:imax+1,0.0)
    slug_temp = splice!(slug_flags,1:imax+1,edge_status)
    append!(Xslugs,[xtemp;_tube_length(sys.tube)])
    append!(slug_flags,[slug_temp;edge_status])

    Interpolations.deduplicate_knots!(Xslugs,move_knots = true)

    sluginterp = linear_interpolation(Xslugs, slug_flags);
end

# overloads to handle the different PHP systems

tube_length(sys::UniformPhpSystem) = _tube_length(sys.tube)
tube_length(sys::VariableAreaPhpSystem) = _tube_length(sys.tube)
tube_length(sys::PhpSystem_nomapping) = _tube_length(sys.tube)

_film_heat_transfer(sys, film_area, positions) =
    Hfilm.(_film_thickness(sys, film_area, positions), Ref(sys.vapor))

_film_thickness(sys, film_area, positions) =
    getδFromδarea.(Ac(sys, positions), Dh(sys, positions), film_area)

Dh(sys::AbstractPHP, positions) = Dh(sys.tube, positions)
Pc(sys::AbstractPHP, positions) = Pc(sys.tube, positions)
Ac(sys::AbstractPHP, positions) = Ac(sys.tube, positions)
