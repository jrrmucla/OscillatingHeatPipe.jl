using Test
using OscillatingHeatPipe
using JLD2

@testset "Results IO" begin
    mktempdir() do dir
        path = joinpath(dir, "history.jld2")
        history = (tube_hist_t=[0.0, 0.1],
                   tube_hist_u=[[1.0], [2.0, 3.0]],
                   tube_hist_θwall=[[290.0], [291.0]],
                   plate_T_hist=[fill(290.0, 2, 3), fill(291.0, 2, 3)])
        @test save_results(path, history) == path
        @test load_results(path) == history
        checkpoint = SimulationResult(Any[history.tube_hist_t...],
            Any[history.tube_hist_u...], Any[history.tube_hist_θwall...],
            Any[], Any[history.plate_T_hist...], nothing, nothing, nothing, nothing)
        checkpoint_path = joinpath(dir, "checkpoint.jld2")
        save_results(checkpoint_path, checkpoint; kind=:checkpoint)
        @test load_results(checkpoint_path) isa SimulationResult
        @test load_results(checkpoint_path; kind=:checkpoint).tube_hist_u == history.tube_hist_u
        @test load_results(checkpoint_path; kind=:history) == history
        @test load_results(path; kind=:history) == history
        @test_throws ArgumentError load_results(path; kind=:checkpoint)
        @test_throws ArgumentError save_results(path, history; kind=:unknown)
        @test_throws ArgumentError save_results(path, history; kind=:checkpoint)
        @test_throws ArgumentError save_results(path, merge(history, (tube_hist_u=[identity],)))
        # Invalid input is rejected before truncating a valid file.
        @test load_results(path) == history
        jldopen(path, "a+") do file
            delete!(file, "format_version")
            file["format_version"] = 999
        end
        @test_throws ArgumentError load_results(path)
        legacy = joinpath(dir, "legacy.jld2")
        jldsave(legacy; SimulationResult=history)
        @test_throws ArgumentError load_results(legacy)
    end
end
