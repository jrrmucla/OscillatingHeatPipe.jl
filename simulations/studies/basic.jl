##
include(joinpath(@__DIR__, "..", "experiments", "baseline.jl"));
include(joinpath(@__DIR__, "write_animations.jl"));

using OscillatingHeatPipe # our main package
using Plots               # for plotting
using JLD2
using Dates
##
title = Dates.format(now(), "yyyymmdd_HHMMSS");
title = "baseline_" * title;
results_dir = "simulations/results";

##
read_path = run_baseline_experiment(title, Δx=0.0006, Δt=5e-4, tf=10.0, fluid="Butane");

## Read Results
Result = load(read_path, "SimulationResult")
write_animations(Result, title, results_dir, 273.0+18.0, 273.0+85.0)
