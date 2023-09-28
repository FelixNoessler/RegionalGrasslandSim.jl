function grazing(sim;
    grazing_half_factor = 1500,
    path = nothing)
    nspecies = 3
    nbiomass = 500

    LD = 2u"ha ^ -1"
    biomass_vec = LinRange(0, 1000, nbiomass)u"kg / ha"
    ρ = [0.9, 1.0, 1.2]

    calc = (;
        biomass_ρ = fill(0.0, nspecies)u"kg / ha",
        grazed_share = fill(0.0, nspecies),
        defoliation = fill(0.0, nspecies)u"kg / (ha * d)")

    grazing_mat = Array{Float64}(undef, nspecies, nbiomass)

    for (i, biomass) in enumerate(biomass_vec)
        calc.defoliation .= 0.0u"kg / (ha * d)"
        sim.Growth.grazing!(;
            calc,
            LD,
            biomass = repeat([biomass], 3),
            ρ,
            grazing_half_factor)
        grazing_mat[:, i] = ustrip.(calc.defoliation)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1],
        xlabel = "Total biomass [green dry mass kg ha⁻¹]",
        ylabel = "Grazed biomass (graz)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title = "")

    for i in 1:nspecies
        lines!(ustrip.(biomass_vec) .* 3, grazing_mat[i, :];
            color = i,
            colorrange = (1, nspecies),
            linewidth = 3, label = "ρ=$(ρ[i])")
    end

    lines!(ustrip.(biomass_vec) .* 3, vec(sum(grazing_mat, dims = 1));
        linewidth = 3,
        color = :grey,
        markersize = 10,
        label = "total")

    axislegend(; framevisible = false, position = :lt)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function grazing_half_factor(; path = nothing)
    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1],
        xlabel = "Total biomass [green dry mass kg ha⁻¹]",
        ylabel = "Grazed biomass (totgraz)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title = "")

    for grazing_half_factor in [750, 1500, 2000]
        x = 0:3000

        LD = 2
        κ = 22
        k_exp = 2
        μₘₐₓ = κ * LD
        h = 1 / μₘₐₓ
        a = 1 / (grazing_half_factor^k_exp * h)
        y = @. a * x^k_exp / (1^k_exp + a * h * x^k_exp)

        lines!(x, y, label = "$grazing_half_factor",
            linewidth = 3,
            color = grazing_half_factor,
            colorrange = (500, 2000))
    end

    axislegend("grazing_half_factor"; framevisible = true, position = :rb)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end
end

function trampling(sim;
    nspecies = 5,
    nLD = 500,
    biomass = fill(100.0, nspecies)u"kg / ha",
    LDs = LinRange(0.0, 4.0, nLD)u"ha^-1",
    trampling_factor = 100,
    path = nothing)
    height = reverse([0.1, 0.2, 0.5, 0.8, 1.0]u"m")

    calc = (;
        trampling_ω = fill(0.0, nspecies)u"ha^-1",
        trampled_biomass = fill(0.0, nspecies)u"kg / ha",
        trampling_high_LD = fill(false, nspecies),
        defoliation = fill(0.0, nspecies)u"kg / (ha * d)")

    trampling_mat_height = Array{Float64}(undef, nspecies, nLD)

    for (i, LD) in enumerate(LDs)
        calc.defoliation .= 0.0u"kg / (ha * d)"
        sim.Growth.trampling!(;
            calc,
            LD,
            biomass,
            height,
            trampling_factor)

        trampling_mat_height[:, i] = ustrip.(calc.defoliation)
    end
    trampling_mat_height = trampling_mat_height ./ 100.0

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1],
        ylabel = "Proportion of biomass that is\nremoved by trampling [d⁻¹]",
        xlabel = "Livestock density [ha⁻¹]",
        title = "Influence of the plant height")
    for i in 1:nspecies
        lines!(ustrip.(LDs), trampling_mat_height[i, :];
            linewidth = 3, label = "height=$(height[i])",
            colormap = :viridis,
            colorrange = (1, nspecies),
            color = i)
    end
    axislegend(; framevisible = false, position = :lt)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function mowing(sim;
    nspecies = 3,
    nbiomass = 3,
    biomass_vec = LinRange(0, 1000, nbiomass)u"kg / ha",
    height = [0.5, 0.3, 0.1]u"m",
    mowing_height = 0.07u"m",
    mowing_mid_days = 30,
    days_since_last_mowing = 100,
    path = nothing)
    calc = (;
        mown_height = fill(0.0, nspecies)u"m",
        mowing_λ = fill(0.0, nspecies),
        defoliation = fill(0.0, nspecies)u"kg / (ha * d)")

    mowing_mat = Array{Float64}(undef, nspecies, nbiomass)

    for (i, biomass) in enumerate(biomass_vec)
        calc.defoliation .= 0.0u"kg / (ha * d)"
        sim.Growth.mowing!(;
            calc,
            mowing_height,
            days_since_last_mowing,
            height,
            biomass,
            mowing_mid_days)

        mowing_mat[:, i] = ustrip.(calc.defoliation)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1],
        xlabel = "Total biomass [green dry mass kg ha⁻¹]",
        ylabel = "Maximal amount of biomass that is\nremoved by mowing (mow)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title = "")

    for i in 1:nspecies
        lines!(ustrip.(biomass_vec) .* nspecies, mowing_mat[i, :];
            linewidth = 3, label = "height=$(height[i])",
            color = i,
            colorrange = (1, nspecies))
    end

    if nspecies <= 5
        lines!(ustrip.(biomass_vec) .* nspecies, vec(sum(mowing_mat, dims = 1));
            color = :grey,
            markersize = 10,
            label = "total")
    end

    axislegend(; framevisible = false, position = :lt)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function mow_factor(;
    path = nothing)
    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1],
        xlabel = "Time since last mowing event [day]\n(days_since_last_mowing)",
        ylabel = "Regrowth of plants (mow_factor)",
        title = "")

    for mowing_mid_days in [20, 40, 100]
        x = 0:200
        y = @. 1 / (1 + exp(-0.05 * (x - mowing_mid_days)))

        lines!(x, y, label = "$mowing_mid_days",
            linewidth = 3,
            color = mowing_mid_days,
            colorrange = (0, 100))
    end

    axislegend("mowing_mid_days"; framevisible = true, position = :rb)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end
