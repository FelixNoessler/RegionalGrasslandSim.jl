function grazing(sim;
    grazing_half_factor=1500,
    path=nothing)

    fig = Figure(; resolution=(700, 400))
    Axis(fig[1,1],
        xlabel="Total biomass [green dry mass kg ha⁻¹]",
        ylabel="Grazed biomass (graz)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title="")

    nspecies = 3
    nbiomass = 500

    LD = 2u"ha ^ -1"
    biomass_vec = LinRange(0, 1000, nbiomass)u"kg / ha"
    ρ = [0.9, 1.0, 1.2]

    grazing_mat = Array{Float64}(undef, nspecies, nbiomass)

    for (i,biomass) in enumerate(biomass_vec)
        graz = sim.Growth.grazing(;
                LD,
                biomass=repeat([biomass], 3),
                ρ,
                grazing_half_factor,
                nspecies
        )
        grazing_mat[:, i] = ustrip.(graz)
    end


    for i in 1:nspecies
        lines!(ustrip.(biomass_vec) .* 3, grazing_mat[i, :];
            color=i,
            colorrange=(1, nspecies),
            linewidth=3, label="ρ=$(ρ[i])")
    end

    lines!(ustrip.(biomass_vec) .* 3, vec(sum(grazing_mat, dims=1));
        linewidth=3,
        color=:grey,
        markersize=10,
        label="total")

    axislegend(; framevisible=false, position=:lt)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end

function grazing_half_factor(; path=nothing)
    fig = Figure(; resolution=(700, 400))
    Axis(fig[1,1],
        xlabel="Total biomass [green dry mass kg ha⁻¹]",
        ylabel="Grazed biomass (totgraz)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title="")

    for grazing_half_factor in [750, 1500, 2000]
        x = 0:3000

        LD = 2
        κ = 22
        k_exp = 2
        μₘₐₓ = κ * LD
        h = 1 / μₘₐₓ
        a = 1 / (grazing_half_factor^k_exp * h)
        y = @. a * x^k_exp  / (1 ^ k_exp + a*h*x^k_exp)

        lines!(x,y,label="$grazing_half_factor",
            linewidth=3,
            color=grazing_half_factor,
            colorrange=(500, 2000))
    end

    axislegend("grazing_half_factor"; framevisible=true, position=:rb)

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
    trampling_factor=100,
    path=nothing)

    fig = Figure(; resolution=(700, 400))

    ## ----------- Influence of CH
    ax1 = Axis(fig[1,1],
        ylabel="Proportion of biomass that is\nremoved by trampling [d⁻¹]",
        xlabel="Livestock density [ha⁻¹]",
        title="Influence of the plant height")
    CH = reverse([0.1, 0.2, 0.5, 0.8, 1.0]u"m")
    trampling_mat_CH = Array{Float64}(undef, nspecies, nLD)
    for (i,LD) in enumerate(LDs)
        trampled_biomass = sim.
            Growth.trampling(; LD, biomass, CH, nspecies, trampling_factor)
            trampling_mat_CH[:, i] = ustrip.(trampled_biomass)
    end
    trampling_mat_CH = trampling_mat_CH ./ 100.0
    for i in 1:nspecies
        lines!(ustrip.(LDs), trampling_mat_CH[i, :];
            linewidth=3, label="CH=$(CH[i])",
            colormap=:viridis,
            colorrange=(1, nspecies),
            color=i)
    end
    axislegend(; framevisible=false, position=:lt)

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
    biomass_vec = LinRange(0, 1000, nbiomass),
    CH = [0.5, 0.3, 0.1]u"m",
    mowing_height=7,
    mowing_mid_days=30,
    days_since_last_mowing=100,
    path=nothing)

    fig = Figure(; resolution=(700, 400))
    Axis(fig[1,1],
        xlabel="Total biomass [green dry mass kg ha⁻¹]",
        ylabel="Maximal amount of biomass that is\nremoved by mowing (mow)\n[green dry mass kg ha⁻¹ d⁻¹]",
        title="")
    mowing_mat = Array{Float64}(undef, nspecies, nbiomass)

    for (i,biomass) in enumerate(biomass_vec)
        mow =
            sim.Growth.mowing(;
                biomass=repeat([biomass], nspecies),
                CH,
                mowing_height,
                mowing_mid_days,
                days_since_last_mowing
        )
        mowing_mat[:, i] = ustrip.(mow)
    end


    for i in 1:nspecies
        lines!(biomass_vec .* nspecies, mowing_mat[i, :];
            linewidth=3, label="CH=$(CH[i])",
            color=i,
            colorrange=(1, nspecies))
    end

    if nspecies <= 5
        lines!(biomass_vec .* nspecies, vec(sum(mowing_mat, dims=1));
            color=:grey,
            markersize=10,
            label="total")
    end

    axislegend(; framevisible=false, position=:lt)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end


function mow_factor(;
    path=nothing)

    fig = Figure(; resolution=(700, 400))
    Axis(fig[1,1],
        xlabel="Time since last mowing event [day]\n(days_since_last_mowing)",
        ylabel="Regrowth of plants (mow_factor)",
        title="")

    for mowing_mid_days in [20, 40, 100]
        x = 0:200
        y = @. 1/(1+exp(-0.05*(x-mowing_mid_days)))

        lines!(x,y,label="$mowing_mid_days",
            linewidth=3,
            color=mowing_mid_days,
            colorrange=(0, 100))
    end

    axislegend("mowing_mid_days"; framevisible=true, position=:rb)

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end
