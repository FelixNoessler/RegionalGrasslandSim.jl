function temperatur_reducer(sim;
    Ts = LinRange(0.0, 40.0, 500), # °C
    path = nothing)
    Ts = sort(ustrip.(Ts))

    y = Float64[]
    for T in Ts
        g = sim.Growth.temperature_reduction(; T)
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
        g = sim.Growth.radiation_reduction(; PAR)
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


function ch_influence(sim;
        path = nothing,
        CH_strength)
    fig = Figure(; resolution = (700, 400))
    ax = Axis(fig[1, 1];
        ylabel = "Plant height growth factor (CHinfluence)",
        xlabel = "Plant height (CH) [m]")

    biomass = fill(10, 3)u"kg / ha"
    CH = [NaN, 0.1, 0.5]u"m"

    CH_vals = 0.0:0.01:1.0
    growth_factors = Array{Float64}(undef, 3, length(CH_vals))



    for (i,CH_1) in enumerate(CH_vals)
        CH[1] = CH_1 * u"m"
        fs = sim.Growth.height_influence(; biomass, CH, CH_strength)
        growth_factors[:, i] .= fs
    end
    ax.title = "CH_strength = $(CH_strength)"
    lines!(CH_vals, growth_factors[1, :];
        label = "varied CH on x-Axis",
        linewidth = 3)
    lines!(CH_vals, growth_factors[2, :];
        label = "CH=$(CH[2])",
        linewidth = 3)
    lines!(CH_vals, growth_factors[3, :];
        label = "CH=$(CH[3])",
        linewidth = 3)
    Legend(fig[1,2], ax; framevisible = false)

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
        g = sim.Growth.seasonal_reduction(; ST)
        push!(y, g)
    end

    fig = Figure(; resolution = (700, 400))
    Axis(fig[1, 1];
        ylabel = "Seasonal factor (SEA)",
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
