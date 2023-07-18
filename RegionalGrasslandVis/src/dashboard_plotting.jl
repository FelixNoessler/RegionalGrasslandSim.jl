
function band_patch(ax; patch=1, sol, color, colormap, colorrange)

    tsub = Int.(round.(LinRange(1, length(sol.t), 1000)))

    ax.ylabel="Green biomass [t ha⁻¹]"
    ax.xlabel="Time [years]"

    biomass_sum = zeros(size(sol.biomass, 1))

    is = sortperm(color)

    cmap = cgrad(colormap)
    colors =[
        cmap[(co .- colorrange[1]) ./ (colorrange[2] - colorrange[1])]
        for co in color[is]
    ]

    for i in eachindex(is)
        biomass_i = ustrip.(uconvert.(u"Mg / ha", sol.biomass[:, patch, is[i]]))

        band!(ax, sol.t[tsub] ./ 365,
            biomass_sum[tsub], (biomass_sum + biomass_i)[tsub];
            color=colors[i]
        )

        biomass_sum += biomass_i
    end

    return nothing
end

function growth_rates(ax; patch=1, sol, color, colormap, colorrange)

    biomass = sol.biomass
    bio1 = biomass[1:end-1, patch, :]
    bio2 = biomass[2:end, patch, :]
    growth = ustrip.((bio2 .- bio1) ./ ((bio2 .+ bio1) ./ 2)  )

    ax.xlabel= "Time [months]"
    ax.ylabel="Net growth [kg/kg ha⁻¹ d⁻¹]"

    skip = length(sol.t) > 1000 ? 10 : 1

    for i in 1:sol.p.nspecies
        lines!(ax, (0.5 .+ sol.t[2:skip:end]) ./ 365 .* 12, growth[1:skip:end, i];
            color=color[i],
            colormap=(colormap, 0.8), colorrange)
    end

    # ylims!(ax, -0.3*maximum(growth), nothing)

    return nothing
end

function trait_mean_biomass(trait, trait_name, ax; sol, color, colormap, colorrange, markersize=15)
    ax.ylabel="Mean biomass [t ha⁻¹]"
    ax.xlabel=trait_name

    scatter!(ax,
        ustrip.(sol.p.species[trait]),
        ustrip.(uconvert.(u"Mg / ha", species_biomass(sol.biomass)));
        color, colormap, colorrange,
        markersize)

    return nothing
end
