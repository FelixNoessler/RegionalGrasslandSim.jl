function optimize_biomass(sim::Module;
    inf_p,
    plotID,
    use_satbiomass = "only",
    startyear = 2012,
    endyear = 2021)
    inf_p = (; inf_p...)
    data, sol = get_plottingdata(sim;
        inf_p,
        plotID,
        startyear,
        endyear)

    diff = 0.0

    if use_satbiomass == "only"
        ### select the days where we have biomass measured
        sim_biomass = @view sol.biomass[data.biomass_t, :, :]

        ### calculate the sum of biomass of all plant species
        biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

        diff = sum(abs.(biomass_sum .- data.biomass))
    elseif use_satbiomass == "no"
        ### select the days where we have biomass measured
        f = 5 * 365 .< data.measured_biomass_t .< 10 * 365
        sim_biomass = @view sol.biomass[data.measured_biomass_t[f], :, :]

        ### calculate the sum of biomass of all plant species
        biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

        diff = sum(abs.(biomass_sum .- data.measured_biomass[f]))
    else
        ### select the days where we have biomass measured
        sim_biomass = @view sol.biomass[data.biomass_t, :, :]

        ### calculate the sum of biomass of all plant species
        biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

        diff += 0.5 * sum(abs.(biomass_sum .- data.biomass))

        ### select the days where we have biomass measured
        f = 5 * 365 .< data.measured_biomass_t .< 10 * 365
        sim_biomass = @view sol.biomass[data.measured_biomass_t[f], :, :]

        ### calculate the sum of biomass of all plant species
        biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

        diff += 0.5 * sum(abs.(biomass_sum .- data.measured_biomass[f]))
    end

    return diff
end
