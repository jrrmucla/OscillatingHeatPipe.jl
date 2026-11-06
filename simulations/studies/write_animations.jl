export write_animations

using OscillatingHeatPipe
using Plots              
using JLD2

function write_animations(SimuResult, title, results_dir, Tref, Tmax)
    t = SimuResult.tube_hist_t

    # Slugs
    anim = @animate for i in eachindex(t)
        plot(OHPSlug(),i,SimuResult)
    end
    gif(anim, joinpath(results_dir, "$(title)_slug.gif"))

    # Temperature
    anim = @animate for i in eachindex(t)
        plot(OHPTemp(),i,SimuResult,clim=(Tref,Tmax))
    end
    gif(anim, joinpath(results_dir, "$(title)_temp.gif"))
end