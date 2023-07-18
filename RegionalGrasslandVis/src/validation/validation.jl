function add_nans(t, b)
    append!(t, [365*i for i in 1:4])
    append!(b, [NaN for i in 1:4])

    return sort(t), b[sortperm(t)]
end


function biomass_validation(data, sol;
        inityears=5, plotID="")

    biomass_sim = vec(sum(ustrip.(sol.biomass); dims=3))

    data_t_prep, data_biomass = add_nans(data.biomass_t, data.biomass)
    data_t = 2012 .+ (inityears * 365 .+ data_t_prep) ./ 365

    mdata_t = 2012 .+ (inityears * 365 .+ data.measured_biomass_t) ./ 365
    mdata = data.measured_biomass
    mdata1 = data.measured_biomass1

    sim_t = 2012 .+ sol.t ./ 365

    fig = Figure()
    Axis(fig[1,1];
        title=plotID,
        xticks=2012:2022,
        xlabel="Time [years]",
        ylabel="Total biomass [kg/ha]",
        limits=(2011, 2023, 0, 5000))
    scatterlines!(
        data_t,
        data_biomass;
        label="Statistical satellite model",
        color=(:red, 0.6),
        linewidth=5)
    lines!(sim_t, biomass_sim;
        label="Mechanistic simulation model",
        color=(:green, 0.6),
        linewidth=5)
    scatter!(mdata_t, mdata,
        marker=:hexagon, color=:black, markersize=15,
        label="Cutted and weighted biomass")
    scatter!(mdata_t, mdata1, marker=:hexagon, color=:black, markersize=15)

    axislegend(; position=:lt)

    return fig
end


function soilmoisture_validation(
        data, sol;
        inityears=5, converted=true, patch=1)

    moisture = ustrip.(sol.water[:, patch]) ./ sol.p.site.root_depth
    moisture[moisture .> 1] .= 1
    moisture[1] = NaN

    if converted
        moisture .*= sol.p.inf_p.moisture_conv
    end

    data_t = 2012 .+ (inityears * 365 .+ (1:length(data.soil_moisture))) ./ 365
    sim_t = 2012 .+ sol.t ./ 365
    fig = Figure()
    Axis(fig[1,1];
        xticks=2012:2022)
    scatterlines!(
        data_t,
        data.soil_moisture;
        label="Measured soil moisture")
    scatterlines!(sim_t, moisture;
        label="Simulated soil moisture$(converted ? "_conv" : "")")

    axislegend(; position=:lb)

    return fig
end


function evaporation_validation(data, sol; inityears=5, patch=1)
    evaporation = ustrip.(sol.evaporation[:, patch])

    data_evaporation = data.evaporation
    data_evaporation[ismissing.(data_evaporation)] .= NaN

    data_t = 2012 .+ (inityears * 365 .+ (1:length(data.evaporation))) ./ 365
    sim_t = 2012 .+ sol.t ./ 365
    fig = Figure()
    Axis(fig[1,1];
        xticks=2012:2022,
        ylabel="Evaporation [mm]",
        xlabel="Time [years]")
    scatterlines!(
        data_t,
        data_evaporation;
        label="Measured evaporation")
    scatterlines!(sim_t, evaporation;
        label="Simulated evaporation")

    axislegend(; position=:lt)

    return fig
end
