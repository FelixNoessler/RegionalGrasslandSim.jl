function amc_nut_response(sim;
    nspecies=6,
    mycorrhizal_colon = LinRange(0, 0.7, nspecies),
    max_AMC_nut_reduction,
    path=nothing
    )

    mycorrhizal_colon = sort(ustrip.(mycorrhizal_colon))

    Ks, x0s = sim.FunctionalResponse.amc_nut_response(;
        mycorrhizal_colon)

    x = 0:0.01:1

    fig = Figure(resolution=(900, 500))
    Axis(fig[1:2, 1];
        xlabel="Nutrient index",
        ylabel="Growth reduction factor\n← no growth, less reduction →",
        title="Influence of the mycorrhizal colonisation")

    for i in eachindex(Ks)
        fun_response = (;
            myco_nutr_right_bound=Ks[i],
            myco_nutr_midpoint=x0s[i])
        y = sim.Growth.amc_nut_reduction(;
            fun_response,
            x,
            max_AMC_nut_reduction)
        lines!(x, y;
            color=i,
            colorrange=(1, nspecies))

        ##### right upper bound
        scatter!([1], [Ks[i]];
            marker=:ltriangle,
            color=i,
            colorrange=(1, nspecies))

        ##### midpoint
        x0_y = Ks[i] / 2
        scatter!([x0s[i]], [x0_y];
            marker=:x,
            color=i,
            colorrange=(1, nspecies)
            )
    end
    ylims!(-0.05, 1.05)

    Axis(fig[1, 2];
        ylabel="Right upper bound",
        xticklabelsvisible=false)
    for i in eachindex(Ks)
        scatter!(mycorrhizal_colon[i], Ks[i];
            marker=:ltriangle,
            color=i,
            colorrange=(1, nspecies))
    end

    Axis(fig[2, 2];
        xlabel="Mycorrhizal colonisation",
        ylabel="Nutrient index\nat midpoint")
    for i in eachindex(x0s)
        scatter!(mycorrhizal_colon[i], x0s[i];
            marker=:x,
            color=i,
            colorrange=(1, nspecies))
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
    SRSA_above = LinRange(0.0, 0.4, nspecies),
    max_SRSA_water_reduction,
    path=nothing)

    SRSA_above = sort(ustrip.(SRSA_above))

    Ks, x0s = sim.
        FunctionalResponse.
        srsa_response(; SRSA_above)

    x = 0:0.01:1

    fig = Figure(resolution=(900, 500))
    Axis(fig[1:2, 1],
        xlabel="Scaled water availability",
        ylabel="Growth reduction factor\n← no growth, less reduction →")

    for (i, (K, x0)) in enumerate(zip(Ks, x0s))
        fun_response = (; srsa_right_bound=K, srsa_midpoint=x0)
        y = sim.Growth.srsa_water_reduction(;
            fun_response, x, max_SRSA_water_reduction)

        lines!(x, y;
            color=i,
            colorrange=(1, nspecies))

        ##### right upper bound
        scatter!([1], [K];
            marker=:ltriangle,
            color=i,
            colorrange=(1, nspecies))

        ##### midpoint
        x0_y = K / 2
        scatter!([x0], [x0_y];
            marker=:x,
            color=i,
            colorrange=(1, nspecies))
    end
    ylims!(-0.1, 1.1)

    Axis(fig[1, 2];
        xticklabelsvisible=false,
        ylabel="Right upper bound")
    scatter!(SRSA_above, Ks;
        marker=:ltriangle,
        color=1:nspecies,
        colorrange=(1, nspecies))

    Axis(fig[2, 2];
        xlabel="SRSA / above ground biomass",
        ylabel="Scaled water availability\nat midpoint")
    scatter!(SRSA_above, x0s;
        marker=:x,
        color=1:nspecies,
        colorrange=(1, nspecies))

    Label(fig[0, 1:2], "Influence of the root surface area / above ground biomass";
        halign=:left,
        font=:bold)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end


function srsa_nut_response(sim;
    nspecies = 6,
    SRSA_above = LinRange(0.0, 0.4, nspecies),
    max_SRSA_nut_reduction,
    path=nothing)

    SRSA_above = sort(ustrip.(SRSA_above))

    Ks, x0s = sim.
        FunctionalResponse.
        srsa_response(; SRSA_above)

    x = 0:0.01:1

    fig = Figure(resolution=(900, 500))
    Axis(fig[1:2, 1],
        xlabel="Nutrient index",
        ylabel="Growth reduction factor\n← no growth, less reduction →")

    for (i, (K, x0)) in enumerate(zip(Ks, x0s))
        fun_response = (; srsa_right_bound=K, srsa_midpoint=x0)
        y = sim.Growth.srsa_nut_reduction(;
            fun_response, x, max_SRSA_nut_reduction)

        lines!(x, y;
            color=i,
            colorrange=(1, nspecies))

        ##### right upper bound
        scatter!([1], [K];
            marker=:ltriangle,
            color=i,
            colorrange=(1, nspecies))

        ##### midpoint
        x0_y = K / 2
        scatter!([x0], [x0_y];
            marker=:x,
            color=i,
            colorrange=(1, nspecies))
    end
    ylims!(-0.1, 1.1)

    Axis(fig[1, 2];
        xticklabelsvisible=false,
        ylabel="Right upper bound")
    scatter!(SRSA_above, Ks;
        marker=:ltriangle,
        color=1:nspecies,
        colorrange=(1, nspecies))

    Axis(fig[2, 2];
        xlabel="SRSA / above ground biomass",
        ylabel="Nutrient index\nat midpoint")
    scatter!(SRSA_above, x0s;
        marker=:x,
        color=1:nspecies,
        colorrange=(1, nspecies))

    Label(fig[0, 1:2], "Influence of the root surface area / above ground biomass";
        halign=:left,
        font=:bold)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end


function potential_growth(sim;
    nspecies=4,
    SLA=LinRange(0.01, 0.04, nspecies)u"m^2/g",
    biomass = repeat([1], nspecies)u"kg / ha",
    PARs = LinRange(0, 12, 10)u"MJ * d^-1 * m^-2",
    path=nothing
    )

    SLA = round.(sort(ustrip.(SLA)); digits=4) .* unit(SLA[1])

    pot_growth = Array{Float64}(undef, nspecies, 10)

    for (i,PAR) in enumerate(PARs)
        growth_par, _ = sim.Growth.potential_growth.(;
            SLA, nspecies, biomass, PAR)

        pot_growth[:, i] .= ustrip.(growth_par)
    end

    fig = Figure(; resolution=(800, 400))
    Axis(fig[1,1],
        xlabel="Photosynthetically active radiation [MJ m⁻² d⁻¹]",
        ylabel="Potential growth per biomass\n[green dry mass kg ha⁻¹ d⁻¹]",
        title="Influence of the specific leaf area (SLA)")

    colormap = :viridis
    sla = ustrip.(SLA)
    colorrange = (minimum(sla), maximum(sla))

    for i in nspecies:-1:1
        lines!(ustrip.(PARs), pot_growth[i, :];
            label="$(sla[i])",
            colormap,
            colorrange,
            color=sla[i])
    end
    axislegend("SLA [m² g⁻¹]"; position=:lt, framevisible=false)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end


function sla_water_response(sim;
    nspecies = 6,
    SLA = LinRange(0.009, 0.05, nspecies),
    max_SLA_water_reduction,
    path=nothing)

    SLA = sort(ustrip.(SLA))

    x0s = sim.
        FunctionalResponse.
        sla_water_response(; SLA)

    x = 0:0.01:1

    fig = Figure(resolution=(900,400))
    Axis(fig[1, 1];
        xlabel="Scaled water availability",
        ylabel="Growth reduction factor\n← no growth, less reduction →",
        title="Influence of the specific leaf area")

    for (i,x0) in enumerate(x0s)
        fun_response = (; sla_water_midpoint=x0)
        y = sim.Growth.sla_water_reduction(;
            fun_response, x,
            max_SLA_water_reduction)
        lines!(x, y;
            color=i,
            colorrange=(1, nspecies))

        ##### midpoint
        x0_y = 0.5
        scatter!([x0], [x0_y];
            marker=:x,
            color=i,
            colorrange=(1, nspecies))
    end

    ylims!(-0.1, 1.1)
    xlims!(-0.02, 1.02)

    Axis(fig[1, 2];
        xlabel="Specific leaf area [m² g⁻¹]",
        ylabel="Scaled water availability\nat midpoint")
    scatter!(SLA, x0s;
        marker=:x,
        color=1:nspecies,
        colorrange=(1, nspecies))

    if !isnothing(path)
        save(path, fig)
    else
        display(fig)
    end

    return nothing
end
