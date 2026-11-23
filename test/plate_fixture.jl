# new state vector: [ζ_L; ζ_R; Q_B; M; Af_L; Af_R; Lfilm_start; Lfilm_end].
println("Building plate fixture, this may take a few minutes...")

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))

hc = 4500 # condenser coeff
h_adiabatic = 160 # adiabatic coeff
length_ohpchannel_mm = 141 # x range of the ohp in mm

ρₛ = 2730; # material density [kg/m^3] Aluminum3003-O
cₛ  = 8.93e02; # material specific heat [J/kg K]
kₛ  = 1.93e02; # material heat conductivity Aluminum3003-O

αₛ = kₛ/ρₛ/cₛ

dₛ = 1.37e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)

Tref = 291.15;
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)

tube_diameter = 1e-3

Δx = 0.000633 # [m]

Lx = 6*INCHES*1.02; # plate size x [m]
Ly = 1.981*INCHES*1.06; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits

g = PhysicalGrid(1.03 .* xlim,1.1 .* ylim,Δx); # build a gird slightly larger than the plate

Δs = 1.4*cellsize(g)

trim = 0.006
cond_block = 1.1INCHES
xbound = [ -Lx/2,-Lx/2, 
            Lx/2, Lx/2];
ybound = [  Ly/2,-Ly/2, 
           -Ly/2, Ly/2];
body = Polygon(xbound,ybound,Δs)

X = MotionTransform([0,0],0)
joint = Joint(X)
m = RigidBodyMotion(joint,body)
x = zero_motion_state(body,m)
update_body!(body,x,m)

Lx_a = 6*INCHES; # plate size x [m]
Ly_a = 1.981*INCHES; # plate size y [m]
xbound_a = [ -Lx_a/2,-Lx_a/2, 
            Lx_a/2, Lx_a/2];
ybound_a = [  Ly_a/2,-Ly_a/2, 
           -Ly_a/2, Ly_a/2];
body_a = Polygon(xbound_a,ybound_a,Δs)
update_body!(body_a,x,m)

function get_qbplus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbplus = zeros_surface(base_cache)
    return qbplus
end

function get_qbminus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbminus = zeros_surface(base_cache)
    # qbminus .= nrm.u
    return qbminus
end

bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)

function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["areal_heat_capacity"] # [K/s]
end

function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["areaheater_coeff"]
    unsteady_scale = phys_params["areal_heat_capacity"] 

    σ .= h*(T0 - T) / unsteady_scale # [K/s]
end

function adiabaticmodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["adiabatic_coeff"]
    corr = phys_params["areal_heat_capacity"] 

    σ .= h*(T0 - T)/ corr
end

fregion1_h = Rectangle(25e-3,25e-3,1.4*Δx) # H01
tr1_h = RigidTransform((0.0,-0.0),0.0)

heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!);

fregion1_c = Rectangle(7.5e-3,Ly_a/2,1.4*Δx)
tr1_c = RigidTransform((3INCHES-7.5e-3,-0.0),0.0)

fregion2_c = deepcopy(fregion1_c)
tr2_c = RigidTransform((-(3INCHES-7.5e-3),-0.0),0.0)

cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!);
cond2 = AreaForcingModel(fregion2_c,tr2_c,adiabaticmodel!);

ds = 1.5Δx
nturn = 16
width_ohp = 1.821*INCHES - tube_diameter
length_ohp = 5.34*INCHES + (length_ohpchannel_mm *1e-3 - 5.34*INCHES)
gap = 0.11INCHES - 0.5*tube_diameter
pitch = width_ohp/(2*nturn+1)

x0, y0 = Lx_a/2 - 0.34*INCHES - 0.5*tube_diameter + (length_ohpchannel_mm *1e-3 - 5.34*INCHES), width_ohp/2
x, y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,3pi/2)

ohp = BasicBody(x,y) # build a BasicBody based on x,y

tr_ohp = RigidTransform((0.0,0.0),0.0)

function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
    σ .= phys_params["ohp_flux"] ./ phys_params["areal_heat_capacity"] 
end
ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!);

forcing_dict = Dict("heating models" => [heater1,cond1,cond2,ohp_linesource])

tspan = (0.0, 300.0); # start time and end time
dt_record = 0.2   # saving time interval

tstep = 4e-4     # actual time marching step

timestep_fixed(u,sys) = tstep

# Each call owns its parameters; tests override only the values they need.
function _construct_plate_fixture(; phys_params=Dict(), forcing=forcing_dict,
                           grid=g, boundary=body, bc=bcdict, motions=m,
                           timestep_func=timestep_fixed)
    parameters = Dict{String,Any}(
        "diffusivity" => αₛ, "areal_heat_capacity" => ρₛ*cₛ*dₛ,
        "Fourier" => 1.0, "ohp_flux" => [NaN],
        "areaheater_power" => nothing, "areaheater_area" => 50e-3 * 50e-3,
        "areaheater_temp" => 0.0, "areaheater_coeff" => nothing,
        "adiabatic_coeff" => nothing, "background temperature" => Tref)
    merge!(parameters, deepcopy(phys_params))
    prob = NeumannHeatConductionProblem(grid, boundary; scaling=GridScaling,
        phys_params=parameters, bc, motions, forcing, timestep_func)
    return construct_system(prob)
end

function get_powered_plate_fixture(; phys_params=Dict(), kwargs...)
    parameters = Dict{String,Any}(
        "areaheater_power" => 30.0,
        "areaheater_coeff" => 4500.0, 
        "adiabatic_coeff"  => 160.0)
    merge!(parameters, phys_params)
    return _construct_plate_fixture(; phys_params=parameters, kwargs...)
end

function get_adiabatic_plate_fixture(; phys_params=Dict(), kwargs...)
    parameters = Dict{String,Any}(
        "areaheater_power" => 0.0,
        "areaheater_coeff" => 0.0, 
        "adiabatic_coeff"  => 0.0)
    merge!(parameters, phys_params)
    return _construct_plate_fixture(; phys_params=parameters, kwargs...)
end