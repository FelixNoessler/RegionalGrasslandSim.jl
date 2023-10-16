function plot_clonalgrowth(sim; path = nothing)
    nspecies = 1
    npatches = 25
    patch_xdim = 5
    patch_ydim = 5

    neighbours = sim.get_allneighbours(;
        p = (; patch_xdim, patch_ydim, npatches))
    surroundings = vcat.(neighbours, Base.OneTo(npatches))

    p = (; neighbours, surroundings, npatches, patch_xdim, patch_ydim)

    calc = (;
        clonalgrowth = fill(NaN, npatches, nspecies)u"kg / ha",
        relbiomass = fill(NaN, npatches), ##needs to be specified
        biomass_per_patch = fill(NaN, npatches)u"kg / ha",
        u_biomass = fill(0.0001, npatches, nspecies)u"kg / ha",
    )
    calc.u_biomass[13, :] .= 10.0u"kg / ha"
    startcondition = copy(ustrip.(calc.u_biomass[:, 1]))

    sim.Growth.calculate_relbiomass!(; calc, p)
    sim.Growth.clonalgrowth!(; p, calc)

    endcondition = copy(ustrip.(calc.u_biomass[:, 1]))

    xs = [[x for x in 1:p.patch_xdim, _ in 1:p.patch_ydim]...]
    ys = [[y for _ in 1:p.patch_xdim, y in 1:p.patch_ydim]...]

    colorrange = quantile(log10.(startcondition), [0.0, 1.0])

    fig = Figure(; resolution = (600, 300))
    ax1 = Axis(fig[1,1]; title = "startcondition")
    plt = scatter!(xs, ys;
        marker = :rect,
        markersize = 1.5,
        markerspace = :data,
        color = log10.(startcondition),
        colorrange,
        colormap = :viridis)

    ax2 = Axis(fig[1,2]; title = "after clonal growth")
    scatter!(xs, ys;
        marker = :rect,
        markersize = 1.5,
        markerspace = :data,
        color = log10.(endcondition),
        colorrange,
        colormap = :viridis)


    for ax in [ax1, ax2]
        ax.aspect = DataAspect()
        ax.yticks = 1:patch_ydim
        ax.xticks = 1:patch_xdim
        ax.limits = (0, patch_xdim + 1, 0, patch_ydim + 1)
    end

    Colorbar(fig[1, 3], plt, label = "log10 biomass [kg ha⁻¹]")

    if !isnothing(path)
        save(path, fig;)
    else
        display(fig)
    end

    return nothing
end
