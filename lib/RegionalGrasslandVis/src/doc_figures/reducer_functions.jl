function temperatur_reducer(sim;
    Ts = LinRange(0.0, 40.0, 500), # °C
    path = nothing)
    Ts = sort(ustrip.(Ts))

    y = Float64[]
    for T in Ts
        g = sim.Growth.temperature_reduction(; T, temperature_red = true)
        push!(y, g)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1];
        ylabel = "Growth reduction",
        xlabel = "Air temperature [°C]",
        title = "Temperature reducer function")

    if length(y) > 500
        scatter!(Ts, y,
            markersize = 5,
            color = (:coral3, 0.5))
    else
        lines!(Ts, y,
            linewidth = 3,
            color = :coral3)
    end

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function radiation_reducer(sim;
    PARs = LinRange(0.0, 15.0, 1000)u"MJ / (m^2 * d)",
    path = nothing)
    PARs = sort(ustrip.(PARs)) .* unit(PARs[1])

    y = Float64[]

    for PAR in PARs
        g = sim.Growth.radiation_reduction(; PAR, radiation_red = true)
        push!(y, g)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1];
        ylabel = "Growth reduction (Rred)",
        xlabel = "Photosynthetically active radiation (PAR) [MJ m⁻² d⁻¹]",
        title = "Radiation reducer function")

    PARs = ustrip.(PARs)

    if length(y) > 1000
        scatter!(PARs, y,
            markersize = 5,
            color = (:magenta, 0.05))
    else
        lines!(PARs, y,
            linewidth = 3,
            color = :magenta)
    end
    ylims!(-0.05, 1.05)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function height_influence(sim;
    path = nothing,
    height_strength)

    biomass = fill(10, 3)u"kg / ha"
    height = [NaN, 0.1, 0.5]u"m"

    height_vals = 0.0:0.01:1.0
    growth_factors = Array{Float64}(undef, 3, length(height_vals))

    nspecies = length(height)
    calc_var = fill(0.0, nspecies)
    biomass_height = fill(0.0, nspecies)u"kg * m / ha"

    for (i, height_1) in enumerate(height_vals)
        height[1] = height_1 * u"m"
        sim.Growth.height_influence!(;
            biomass,
            height,
            biomass_height,
            height_strength, height_included = true,
            height_influence = calc_var)
        growth_factors[:, i] .= calc_var
    end

    fig = Figure(; resolution = (700, 400))
    ax = Axis(fig[1, 1];
        ylabel = "Plant height growth factor (heightinfluence)",
        xlabel = "Plant height [m]",
        xticklabelcolor = :blue)
    ax.title = "height_strength = $(height_strength)"
    lines!(height_vals, growth_factors[1, :];
        label = "varied height on x-Axis",
        linewidth = 3,
        color = :blue)
    lines!(height_vals, growth_factors[2, :];
        label = "height=$(height[2])",
        linewidth = 3,
        color = :red)
    lines!(height_vals, growth_factors[3, :];
        label = "height=$(height[3])",
        linewidth = 3,
        color = :orange)
    Legend(fig[1, 2], ax; framevisible = false)

    ylims!(ax, 0.1, 2.5)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function below_influence(sim;
    path = nothing,
    below_competition_strength)

    biomass_vals = LinRange(0, 500, 100)u"kg / ha"
    trait_similarity = [1 0.7 0.8;
                        0.7 1 0.75;
                        0.8 0.75 1]u"ha / kg"

    nspecies = size(trait_similarity, 1)
    biomass = fill(25.0, nspecies)u"kg / ha"

    growth_factors = Array{Float64}(undef, nspecies, length(biomass_vals))

    calc_var = fill(0.0, nspecies)
    traitsimilarity_biomass = fill(0.0, nspecies)


    for (i, biomass_val) in enumerate(biomass_vals)
        biomass[1] = biomass_val
        sim.Growth.below_ground_competition!(;
            below = calc_var,
            traitsimilarity_biomass,
            biomass, below_included=true,
            trait_similarity,
            below_competition_strength)
        growth_factors[:, i] .= calc_var
    end

    fig = Figure(; resolution = (700, 400))
    ax = Axis(fig[1, 1];
        ylabel = "Below ground competition\ngrowth factor (below)",
        xlabel = "Biomas of species 1 [kg ha⁻¹]",
        xticklabelcolor = :blue)
    ax.title = "below_competition_strength = $(below_competition_strength)"

    lines!(ustrip.(biomass_vals), growth_factors[1, :];
        label = "varied biomass on x-Axis",
        linewidth = 3,
        color = :blue)
    lines!(ustrip.(biomass_vals), growth_factors[2, :];
        label = "biomass=$(biomass[2])",
        linewidth = 3,
        color = :red)
    lines!(ustrip.(biomass_vals), growth_factors[3, :];
        label = "biomass=$(biomass[3])",
        linewidth = 3,
        color = :orange)
    Legend(fig[1, 2], ax; framevisible = false)

    ylims!(ax, 0.1, 2.5)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function seasonal_effect(sim;
    STs = LinRange(0, 3500, 1000),
    path = nothing)
    y = Float64[]
    for ST in STs
        g = sim.Growth.seasonal_reduction(; ST, season_red = true)
        push!(y, g)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1];
        ylabel = "Seasonal factor (seasonal)",
        xlabel = "Accumulated degree days (ST) [°C]",
        title = "Seasonal effect")

    if length(y) > 1000
        scatter!(STs, y;
            markersize = 3,
            color = (:navajowhite4, 0.1))
    else
        lines!(STs, y;
            linewidth = 3,
            color = :navajowhite4)
    end

    ylims!(-0.05, 1.6)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function seasonal_component_senescence(sim;
    STs = LinRange(0, 3500, 1000),
    path = nothing)
    STs = sort(STs)

    y = Float64[]
    for ST in STs
        g = sim.Growth.seasonal_component_senescence(; ST)
        push!(y, g)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1];
        ylabel = "Seasonal factor (SEN)",
        xlabel = "Accumulated degree days (ST) [°C]",
        title = "")

    if length(y) > 1000
        scatter!(STs, y;
            markersize = 3,
            color = (:navajowhite4, 0.1))
    else
        lines!(STs, y;
            linewidth = 3,
            color = :navajowhite4)
    end
    ylims!(-0.05, 3.5)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end
