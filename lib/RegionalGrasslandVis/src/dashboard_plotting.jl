function band_patch(ax;
    patch = 1,
    sol, color, colormap, colorrange, plotID, valid,
    show_bands, show_grazmow, show_validdata)

    # tsub = Int.(round.(LinRange(1, length(sol.t), 5000)))

    ax.ylabel = "Green biomass [kg ha⁻¹]"
    ax.xlabel = "Time [years]"
    is = sortperm(color)

    cmap = cgrad(colormap)
    colors = [cmap[(co .- colorrange[1]) ./ (colorrange[2] - colorrange[1])]
              for co in color[is]]

    t = sol.t ./ 365
    years = 0:(round(sol.t[end] / 365) - 1)

    if !isnothing(plotID)
        t .+= 2012
        ax.xticks = 2012:2:2022
        ax.xminorticks = 2012:1:2022
        ax.xminorticksvisible = true
        years = 2012:2021
    end

    biomass_sum = vec(sum(ustrip.(sol.biomass); dims = 3))
    if show_grazmow
        ymax = maximum(biomass_sum) * 1.5
        add_grazing!(ax, sol; t, ymax)
        add_mowing!(ax, sol; years, ymax)
    end

    if show_bands
        biomass_cumsum = zeros(size(sol.biomass, 1))
        for i in eachindex(is)
            biomass_i = ustrip.(uconvert.(u"kg / ha", sol.biomass[:, patch, is[i]]))
            band!(ax, t,
                biomass_cumsum, (biomass_cumsum + biomass_i);
                color = colors[i])
            biomass_cumsum += biomass_i
        end
    end

    lines!(ax,
        t, biomass_sum;
        color = :orange)

    if !isnothing(plotID) .&& show_validdata
        data = valid.get_validation_data(; plotID)
        mdata = data.measured_biomass
        mdata1 = data.measured_biomass1
        mdata_t = 2012 .+ data.measured_biomass_t ./ 365
        f = mdata_t .< 2022
        scatter!(ax, mdata_t[f], mdata[f],
            marker = :hexagon, color = :black, markersize = 15)
        scatter!(ax, mdata_t[f], mdata1[f],
            marker = :hexagon, color = :black, markersize = 15)

        data_t_prep, data_biomass = add_nans(data.biomass_t, data.biomass)
        data_t = 2012 .+ data_t_prep ./ 365

        scatterlines!(ax, data_t, data_biomass;
            color = :grey)
    end

    return nothing
end

function growth_rates(ax; patch = 1, sol, color, colormap, colorrange, plotID)
    biomass = sol.biomass
    bio1 = biomass[1:(end - 1), patch, :]
    bio2 = biomass[2:end, patch, :]
    growth = ustrip.((bio2 .- bio1) ./ ((bio2 .+ bio1) ./ 2))

    ax.xlabel = "Time [years]"
    ax.ylabel = "Net growth [kg/kg ha⁻¹ d⁻¹]"

    # skip = length(sol.t) > 1000 ? 10 : 1
    skip = 1

    t = (0.5 .+ sol.t[2:skip:end]) ./ 365

    if !isnothing(plotID)
        t .+= 2012
        ax.xticks = 2012:2:2022
        ax.xminorticks = 2012:1:2022
        ax.xminorticksvisible = true
    end

    for i in 1:(sol.p.nspecies)
        lines!(ax, t, growth[1:skip:end, i];
            color = color[i],
            colormap = (colormap, 0.8), colorrange)
    end

    # ylims!(ax, -0.3*maximum(growth), nothing)

    return nothing
end

function trait_mean_biomass(trait,
    trait_name,
    ax;
    sol,
    color,
    colormap,
    colorrange,
    markersize = 15)
    ax.ylabel = "Mean biomass [t ha⁻¹]"
    ax.xlabel = trait_name

    scatter!(ax,
        ustrip.(sol.p.species[trait]),
        ustrip.(uconvert.(u"Mg / ha", species_biomass(sol.biomass)));
        color, colormap, colorrange,
        markersize)

    return nothing
end
