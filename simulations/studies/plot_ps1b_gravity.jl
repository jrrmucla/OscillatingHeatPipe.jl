# Run from the repository root: julia --project=. simulations/studies/plot_ps1b_gravity.jl
include(joinpath(@__DIR__, "..", "experiments", "ps1b.jl"))

# Reproduce the experiment's routing and placement without time marching.
ohp, transform = _ps1b_ohp_geometry(0.0006)
shape = transform(ohp)
@assert all(x -> 0 <= x <= 0.0692, shape.x) "OHP extends outside the plate in x"
@assert all(y -> 0 <= y <= 0.1205, shape.y) "OHP extends outside the plate in y"
potential = ps1b_gravitational_potential(nothing, shape.x, shape.y)
p = plot(OHPGravity(), shape.x, shape.y, potential;
    closedornot=true, title="PS1b gravitational potential", size=(750, 850))
mkpath(DEFAULT_RESULTS_DIR)
output = joinpath(DEFAULT_RESULTS_DIR, "ps1b_gravitational_potential.png")
savefig(p, output)
println(output)

# Extend the prescribed field over the physical plate, without the grid buffer.
plate_x = range(0.0, 0.0692; length=347)
plate_y = range(0.0, 0.1205; length=604)
plate_potential = [ps1b_gravitational_potential(nothing, x, y)
                   for y in plate_y, x in plate_x]
p_plate = plot(OHPGravity(), plate_x, plate_y, plate_potential;
    title="Unfolded Gravitational Potential", size=(750, 850),
    xlims=extrema(plate_x), ylims=extrema(plate_y))
plot!(p_plate, [shape.x; shape.x[1]], [shape.y; shape.y[1]];
    color=:white, linewidth=1.2, label=false)
plate_output = joinpath(DEFAULT_RESULTS_DIR, "ps1b_plate_gravitational_potential.png")
savefig(p_plate, plate_output)
println(plate_output)
