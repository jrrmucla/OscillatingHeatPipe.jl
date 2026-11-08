include(joinpath(@__DIR__, "baseline.jl"))

export plot_alturnating_layout

"""
    plot_alturnating_layout(title="alturnating"; Δx=0.0006)

Copy the baseline plate, heater, condenser, and OHP layout with rectangular
channels of 1 mm depth and alternating 1 mm / 1.5 mm widths. Save a layout with
line width proportional to hydraulic diameter
layout to `simulations/results/<title>_layout.png` and return the plot, tube,
centerline coordinates, and output path. No simulation system is initialized.
"""
function plot_alturnating_layout(title="alturnating"; Δx=0.0006)
    # Match the baseline plate outline and placement of heating/cooling regions.
    plate_length = 0.1524
    plate_width = 0.0648
    condenser_block_length = 1.1INCHES
    trim = 0.006
    left_edge, right_edge = -plate_length / 2, plate_length / 2
    bottom_edge, top_edge = -plate_width / 2, plate_width / 2
    left_shoulder = left_edge + condenser_block_length
    right_shoulder = right_edge - condenser_block_length
    plate_x = [left_edge, left_edge, left_shoulder, left_shoulder,
               right_shoulder, right_shoulder, right_edge, right_edge,
               right_shoulder, right_shoulder, left_shoulder, left_shoulder]
    plate_y = [top_edge, bottom_edge, bottom_edge, bottom_edge + trim,
               bottom_edge + trim, bottom_edge, bottom_edge, top_edge,
               top_edge, top_edge - trim, top_edge - trim, top_edge]
    plate = Polygon(plate_x, plate_y, 1.4Δx)

    physical_parameters = Dict{String,Any}()
    heater_regions, heater_transforms, _ =
        _place_baseline_heaters!(physical_parameters, 40.0, Δx)
    cooler_regions, cooler_transforms, _ =
        _place_baseline_coolers!(physical_parameters, Δx)

    channel_depth = 1e-3
    narrow_width, wide_width = 1e-3, 1.5e-3
    rectangular_area(width) = channel_depth * width
    rectangular_perimeter(width) = 2 * (channel_depth + width)
    hydraulic_diameter(width) = 4 * rectangular_area(width) / rectangular_perimeter(width)

    # Preserve the ASETS centerline used by the baseline experiment.
    turn_count = 16
    ohp_width = 46.25e-3
    straightaway_length = 133.83e-3
    pitch = ohp_width / (2turn_count + 1)
    closure_gap = 1e-3
    x_origin = -straightaway_length / 2 - 2e-3
    y_origin = -ohp_width / 2
    x, y, xf, yf, tube = construct_alternating_ohp_curve(
        turn_count, pitch, straightaway_length, closure_gap, 1.5Δx,
        x_origin, y_origin, false, false, pi / 2,
        hydraulic_diameter(narrow_width), rectangular_perimeter(narrow_width),
        rectangular_area(narrow_width),
        hydraulic_diameter(wide_width), rectangular_perimeter(wide_width),
        rectangular_area(wide_width))

    layout_plot = plot(plate; fillalpha=0, linecolor=:black, legend=false,
                       aspect_ratio=:equal, size=(1000, 500),
                       title="Alternating rectangular channels",
                       xlabel="x [m]", ylabel="y [m]")
    for (region, transform) in zip(heater_regions, heater_transforms)
        update_body!(region, transform)
        plot!(layout_plot, region; fillcolor=:red, fillalpha=0.25, label=false)
    end
    for (region, transform) in zip(cooler_regions, cooler_transforms)
        update_body!(region, transform)
        plot!(layout_plot, region; fillcolor=:blue, fillalpha=0.25, label=false)
    end
    plot!(layout_plot, OHPDh(), x, y, xf, yf, tube)
    plot!(layout_plot; margin=6Plots.mm)

    results_directory = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_directory)
    output_path = joinpath(results_directory, "$(title)_layout.png")
    savefig(layout_plot, output_path)
    return (; layout_plot, tube, x, y, xf, yf, output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    layout = plot_alturnating_layout()
    display(layout.layout_plot)
    println("Saved layout to ", layout.output_path)
end
