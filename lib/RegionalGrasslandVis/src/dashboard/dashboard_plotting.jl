function band_patch(;
    plot_obj,
    patch = 1,
    sol,
    valid_data,
    ax_num = 1)
    ax = plot_obj.axes[ax_num]
    empty!(ax)

    ax.ylabel = "Green biomass [kg ha⁻¹]"
    ax.xlabel = "Time [years]"

    t = sol.numeric_date
    biomass_sum = vec(sum(ustrip.(sol.biomass); dims = 3))

    show_grazmow = plot_obj.obs.toggle_grazmow.active.val
    if show_grazmow
        ymax = maximum(biomass_sum) * 1.5
        add_grazing!(ax, sol; ymax)
        add_mowing!(ax, sol; ymax)
    end

    show_bands = plot_obj.obs.toggle_bands.active.val
    if show_bands
        trait = plot_obj.obs.menu_color.selection.val
        color = ustrip.(sol.p.species[trait])
        colormap = :viridis
        colorrange = (minimum(color), maximum(color) + 0.0001)
        is = sortperm(color)
        cmap = cgrad(colormap)
        colors = [cmap[(co .- colorrange[1]) ./ (colorrange[2] - colorrange[1])]
                  for co in color[is]]

        biomass_cumsum = zeros(size(sol.biomass, 1))
        for i in eachindex(is)
            biomass_i = ustrip.(uconvert.(u"kg / ha", sol.biomass[:, patch, is[i]]))
            band!(ax, t,
                biomass_cumsum, (biomass_cumsum + biomass_i);
                color = colors[i])
            biomass_cumsum += biomass_i
        end

    else
        lines!(ax,
            t, biomass_sum;
            color = :orange)
    end

    if !isnothing(valid_data)
        mdata = values(valid_data.measured_biomass.biomass)
        mdata_t = values(valid_data.measured_biomass.numeric_date)
        scatter!(ax, mdata_t, mdata, color = :black, markersize = 8)
    end

    return nothing
end

function trait_time_plot(;
    sol, patch = 1, valid_data, plot_obj, ax_num = 2)
    ax = plot_obj.axes[ax_num]
    empty!(ax)
    t = sol.numeric_date

    trait = plot_obj.obs.menu_color.selection.val
    name_index = getindex.([plot_obj.obs.menu_color.options.val...], 2) .== trait
    trait_name = first.([plot_obj.obs.menu_color.options.val...])[name_index][1]

    trait_vals = ustrip.(sol.p.species[trait])

    biomass_vals = ustrip.(sol.biomass[:, patch, :])
    total_biomass = sum(biomass_vals, dims = 2)
    relative_biomass = biomass_vals ./ total_biomass

    ##  mean
    weighted_trait = trait_vals .* relative_biomass'
    cwm_trait = vec(sum(weighted_trait; dims = 1))

    ax.xlabel = "Time [years]"

    show_traitvar = plot_obj.obs.toggle_traitvar.active.val
    if show_traitvar
        ## variance
        trait_diff = (trait_vals' .- cwm_trait) .^ 2
        weighted_trait_diff = trait_diff .* relative_biomass
        cwv_trait = vec(sum(weighted_trait_diff; dims = 2))

        lines!(ax, t, cwv_trait, color = :red)
        ax.ylabel = "Var: $trait_name"
    else
        ### trait values of all species
        for i in 1:(sol.p.nspecies)
            trait_i = trait_vals[i]
            lines!(ax, [t[1], t[end]], [trait_i, trait_i], color = (:grey, 0.2))
        end
        lines!(ax, t, cwm_trait, color = :blue)
        ax.ylabel = "Mean: $trait_name"
    end

    if !isnothing(valid_data) && !show_traitvar
        x = valid_data.traits.num_t
        y = vec(valid_data.traits.val[:, trait .== valid_data.traits.dim])
        scatter!(ax, x, y, color = :black, markersize = 8)
    end
end

function trait_mean_biomass(;
    sol,
    markersize = 15,
    patch = 1,
    t, plot_obj, ax_num = 3)
    ax = plot_obj.axes[ax_num]
    empty!(ax)

    trait = plot_obj.obs.menu_color.selection.val
    name_index = getindex.([plot_obj.obs.menu_color.options.val...], 2) .== trait
    trait_name = first.([plot_obj.obs.menu_color.options.val...])[name_index][1]
    color = ustrip.(sol.p.species[trait])
    colormap = :viridis
    colorrange = (minimum(color), maximum(color) + 0.0001)

    print_date = Dates.format(sol.date[t], "dd.mm.yyyy")

    ax.ylabel = "Biomass at $(print_date) [kg ha⁻¹]"
    ax.xlabel = trait_name
    ax.xzoomlock = true
    ax.yzoomlock = true
    ax.xrectzoom = false
    ax.yrectzoom = false
    ax.xpanlock = true
    ax.ypanlock = true

    scatter!(ax,
        ustrip.(sol.p.species[trait]),
        ustrip.(sol.biomass[t, patch, :]);
        color, colormap, colorrange,
        markersize)

    return nothing
end

function soilwater_plot(; sol, valid_data, plot_obj, ax_num = 4)
    ax = plot_obj.axes[ax_num]
    empty!(ax)
    scatterlines!(ax, sol.numeric_date, ustrip.(sol.water[:, 1]);
        color = :turquoise3,
        markersize = 4,
        linewidth = 0.1)
    PWP = ustrip(sol.p.site.PWP)
    WHC = ustrip(sol.p.site.WHC)
    lines!(ax, [sol.numeric_date[1], sol.numeric_date[end]], [PWP, PWP];
        color = :blue)
    lines!(ax, [sol.numeric_date[1], sol.numeric_date[end]], [WHC, WHC];
        color = :blue)
    ax.ylabel = "Soil water [mm]"
    ax.xlabel = "Time [years]"
    ylims!(ax, 0.0, nothing)

    if !isnothing(valid_data)
        α = sol.p.inf_p.moistureconv_alpha
        β = sol.p.inf_p.moistureconv_beta
        scatter!(ax,
            valid_data.soilmoisture.num_t,
            α .+ β .* valid_data.soilmoisture.val;
            color = :black,
            markersize = 2)
    end
end

function abiotic_plot(; sol, plot_obj, ax_num = 6)
    ax = plot_obj.axes[ax_num]
    empty!(ax)
    abiotic_colors = [:blue, :brown, :red, :red, :orange]
    abiotic = plot_obj.obs.menu_abiotic.selection.val
    name_index = getindex.([plot_obj.obs.menu_abiotic.options.val...], 2) .== abiotic
    abiotic_name = first.([plot_obj.obs.menu_abiotic.options.val...])[name_index][1]
    abiotic_color = abiotic_colors[name_index][1]

    scatterlines!(ax, sol.numeric_date, ustrip.(sol.p.daily_data[abiotic]);
        color = abiotic_color,
        markersize = 4,
        linewidth = 0.1)
    ax.ylabel = abiotic_name
    ax.xlabel = "Time [years]"
end

# function growth_rates(ax; patch = 1, sol, color, colormap, colorrange, plotID)
#     biomass = sol.biomass
#     bio1 = biomass[1:(end - 1), patch, :]
#     bio2 = biomass[2:end, patch, :]
#     growth = ustrip.((bio2 .- bio1) ./ ((bio2 .+ bio1) ./ 2))

#     ax.xlabel = "Time [years]"
#     ax.ylabel = "Net growth [kg/kg ha⁻¹ d⁻¹]"

#     # skip = length(sol.t) > 1000 ? 10 : 1
#     skip = 1
#     t = sol.numeric_date[2:end] .+ (0.5 / 365)

#     for i in 1:(sol.p.nspecies)
#         lines!(ax, t, growth[1:skip:end, i];
#             color = color[i],
#             colormap = (colormap, 0.8), colorrange)
#     end

#     # ylims!(ax, -0.3*maximum(growth), nothing)

#     return nothing
# end
