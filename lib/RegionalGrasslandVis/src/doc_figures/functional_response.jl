function amc_nut_response(sim;
    nspecies = 6,
    mycorrhizal_colon = LinRange(0, 0.7, nspecies),
    max_AMC_nut_reduction,
    path = nothing)
    A, Ks, x0s = sim.FunctionalResponse.amc_nut_response(;
        mycorrhizal_colon = sort(mycorrhizal_colon),
        maximal_reduction = max_AMC_nut_reduction)

    calc_var = fill(NaN, nspecies)
    xs = 0:0.01:1
    ymat = fill(0.0, length(xs), nspecies)

    fun_response = (;
        myco_nut_midpoint = x0s,
        myco_nut_lower = A,
        myco_nut_upper = Ks)

    for (i, x) in enumerate(xs)
        sim.Growth.amc_nut_reduction!(;
            amc_nut = calc_var,
            fun_response,
            x)
        ymat[i, :] .= calc_var
    end

    fig = Figure(resolution = (900, 500))
    Axis(fig[1:2, 1];
        xlabel = "Nutrient index",
        ylabel = "Growth reduction factor\n← no growth, less reduction →",
        title = "Influence of the mycorrhizal colonisation")

    for i in eachindex(Ks)
        lines!(xs, ymat[:, i];
            color = i,
            colorrange = (1, nspecies))

        ##### right upper bound
        scatter!([1], [Ks[i]];
            marker = :ltriangle,
            color = i,
            colorrange = (1, nspecies))

        ##### midpoint
        x0_y = (Ks[i] - A) / 2 + A
        scatter!([x0s[i]], [x0_y];
            marker = :x,
            color = i,
            colorrange = (1, nspecies))
    end
    ylims!(-0.05, 1.05)

    Axis(fig[1, 2];
        ylabel = "Right upper bound",
        xticklabelsvisible = false)
    for i in eachindex(Ks)
        scatter!(mycorrhizal_colon[i], Ks[i];
            marker = :ltriangle,
            color = i,
            colorrange = (1, nspecies))
    end

    Axis(fig[2, 2];
        xlabel = "Mycorrhizal colonisation",
        ylabel = "Nutrient index\nat midpoint")
    for i in eachindex(x0s)
        scatter!(mycorrhizal_colon[i], x0s[i];
            marker = :x,
            color = i,
            colorrange = (1, nspecies))
    end

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function srsa_water_response(sim;
    nspecies = 6,
    SRSA_above = LinRange(0.0, 0.4, nspecies)u"m^2 / g",
    max_SRSA_water_reduction,
    path = nothing)
    SRSA_above = collect(sort(SRSA_above))

    A, Ks, x0s = sim.FunctionalResponse.srsa_response(;
        SRSA_above,
        maximal_reduction = max_SRSA_water_reduction)

    calc_var = fill(NaN, nspecies)
    xs = 0:0.01:1
    ymat = fill(0.0, length(xs), nspecies)

    fun_response = (;
        srsa_midpoint = x0s,
        srsa_water_lower = A,
        srsa_water_upper = Ks)

    for (i, x) in enumerate(xs)
        sim.Growth.srsa_water_reduction!(;
            srsa_water = calc_var,
            fun_response,
            x)
        ymat[i, :] .= calc_var
    end

    fig = Figure(resolution = (900, 500))
    Axis(fig[1:2, 1],
        xlabel = "Scaled water availability",
        ylabel = "Growth reduction factor\n← no growth, less reduction →")

    for (i, (K, x0)) in enumerate(zip(Ks, x0s))
        lines!(xs, ymat[:, i];
            color = i,
            colorrange = (1, nspecies))

        ##### right upper bound
        scatter!([1], [K];
            marker = :ltriangle,
            color = i,
            colorrange = (1, nspecies))

        ##### midpoint
        x0_y = (K - A) / 2 + A
        scatter!([x0], [x0_y];
            marker = :x,
            color = i,
            colorrange = (1, nspecies))
    end
    ylims!(-0.1, 1.1)

    Axis(fig[1, 2];
        xticklabelsvisible = false,
        ylabel = "Right upper bound")
    scatter!(ustrip.(SRSA_above), Ks;
        marker = :ltriangle,
        color = 1:nspecies,
        colorrange = (1, nspecies))

    Axis(fig[2, 2];
        xlabel = "SRSA / above ground biomass",
        ylabel = "Scaled water availability\nat midpoint")
    scatter!(ustrip.(SRSA_above), x0s;
        marker = :x,
        color = 1:nspecies,
        colorrange = (1, nspecies))

    Label(fig[0, 1:2], "Influence of the root surface area / above ground biomass";
        halign = :left,
        font = :bold)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function srsa_nut_response(sim;
    nspecies = 6,
    SRSA_above = LinRange(0.0, 0.4, nspecies)u"m^2 / g",
    max_SRSA_nut_reduction,
    path = nothing)
    A, Ks, x0s = sim.FunctionalResponse.srsa_response(;
        SRSA_above = sort(SRSA_above),
        maximal_reduction = max_SRSA_nut_reduction)

    calc_var = fill(NaN, nspecies)
    xs = 0:0.01:1
    ymat = fill(0.0, length(xs), nspecies)

    fun_response = (;
        srsa_midpoint = x0s,
        srsa_nut_lower = A,
        srsa_nut_upper = Ks)

    for (i, x) in enumerate(xs)
        sim.Growth.srsa_nut_reduction!(;
            srsa_nut = calc_var,
            fun_response,
            x)
        ymat[i, :] .= calc_var
    end

    fig = Figure(resolution = (900, 500))
    Axis(fig[1:2, 1],
        xlabel = "Nutrient index",
        ylabel = "Growth reduction factor\n← no growth, less reduction →")

    for (i, (K, x0)) in enumerate(zip(Ks, x0s))
        lines!(xs, ymat[:, i];
            color = i,
            colorrange = (1, nspecies))

        ##### right upper bound
        scatter!([1], [K];
            marker = :ltriangle,
            color = i,
            colorrange = (1, nspecies))

        ##### midpoint
        x0_y = (K - A) / 2 + A
        scatter!([x0], [x0_y];
            marker = :x,
            color = i,
            colorrange = (1, nspecies))
    end
    ylims!(-0.1, 1.1)

    Axis(fig[1, 2];
        xticklabelsvisible = false,
        ylabel = "Right upper bound")
    scatter!(ustrip.(SRSA_above), Ks;
        marker = :ltriangle,
        color = 1:nspecies,
        colorrange = (1, nspecies))

    Axis(fig[2, 2];
        xlabel = "SRSA / above ground biomass",
        ylabel = "Nutrient index\nat midpoint")
    scatter!(ustrip.(SRSA_above), x0s;
        marker = :x,
        color = 1:nspecies,
        colorrange = (1, nspecies))

    Label(fig[0, 1:2], "Influence of the root surface area / above ground biomass";
        halign = :left,
        font = :bold)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function potential_growth(sim;
    nspecies = 4,
    SLA = LinRange(0.01, 0.04, nspecies)u"m^2/g",
    biomass = repeat([1], nspecies)u"kg / ha",
    PARs = LinRange(0, 12, 10)u"MJ * d^-1 * m^-2",
    path = nothing)
    SLA = round.(sort(ustrip.(SLA)); digits = 4) .* unit(SLA[1])
    pot_growth_final = Array{Float64}(undef, nspecies, 10)
    LAIs = Array{Float64}(undef, nspecies)
    pot_growth = fill(NaN, nspecies)u"kg / (ha * d)"

    calc = (; LAIs, pot_growth)

    for (i, PAR) in enumerate(PARs)
        sim.Growth.potential_growth!(;
            calc,
            potgrowth_included = true,
            SLA,
            biomass,
            PAR)

        pot_growth_final[:, i] .= ustrip.(calc.pot_growth)
    end

    fig = Figure(; resolution = (800, 400))
    Axis(fig[1, 1],
        xlabel = "Photosynthetically active radiation [MJ m⁻² d⁻¹]",
        ylabel = "Potential growth per biomass\n[green dry mass kg ha⁻¹ d⁻¹]",
        title = "Influence of the specific leaf area (SLA)")

    colormap = :viridis
    sla = ustrip.(SLA)
    colorrange = (minimum(sla), maximum(sla))

    for i in nspecies:-1:1
        lines!(ustrip.(PARs), pot_growth_final[i, :];
            label = "$(sla[i])",
            colormap,
            colorrange,
            color = sla[i])
    end
    axislegend("SLA [m² g⁻¹]"; position = :lt, framevisible = false)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function sla_water_response(sim;
    nspecies = 6,
    SLA = LinRange(0.009, 0.05, nspecies)u"m^2 / g",
    max_SLA_water_reduction,
    path = nothing)
    A, x0s = sim.FunctionalResponse.sla_water_response(;
        SLA = sort(SLA), maximal_reduction = max_SLA_water_reduction)
    fun_response = (; sla_water_midpoint = x0s, sla_water_lower = A)

    xs = 0:0.01:1
    sla_water = fill(NaN, nspecies)

    fig = Figure(resolution = (900, 400))
    Axis(fig[1, 1];
        xlabel = "Scaled water availability",
        ylabel = "Growth reduction factor\n← no growth, less reduction →",
        title = "Influence of the specific leaf area")

    ymat = fill(0.0, length(xs), nspecies)

    for (i, x) in enumerate(xs)
        sim.Growth.sla_water_reduction!(;
            sla_water,
            fun_response,
            x)
        ymat[i, :] .= sla_water
    end

    for i in eachindex(x0s)
        lines!(xs, ymat[:, i];
            color = i,
            colorrange = (1, nspecies))

        ##### midpoint
        x0_y = 1 - max_SLA_water_reduction / 2
        scatter!([x0s[i]], [x0_y];
            marker = :x,
            color = i,
            colorrange = (1, nspecies))
    end

    ylims!(-0.1, 1.1)
    xlims!(-0.02, 1.02)

    Axis(fig[1, 2];
        xlabel = "Specific leaf area [m² g⁻¹]",
        ylabel = "Scaled water availability\nat midpoint")
    scatter!(ustrip.(SLA), x0s;
        marker = :x,
        color = 1:nspecies,
        colorrange = (1, nspecies))

    if !isnothing(path)
        save(path, fig)
    else
        display(fig)
    end

    return nothing
end
