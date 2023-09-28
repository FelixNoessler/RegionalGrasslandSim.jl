import RegionalGrasslandSim as sim
using RegionalGrasslandValid
using RegionalGrasslandVis ## for plotting theme
using Statistics
using GLMakie
using TimeSeries
using Distributions
Makie.inline!(true)

########################### input preparataion
param_vals = [
    1155.299246041266,
    2.1787199839476785, 11.452466470742415,
    0.381072247902645, 150.1232967984496, 958.214969168516,
    17.821776550046422, 0.9416295869463482, 0.9290589640747696,
    0.6261919010370469, 0.708332191866904]
param_names = [
    "sigma_biomass",
    "senescence_intercept", "senescence_rate",
    "below_competition_strength", "trampling_factor", "grazing_half_factor",
    "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
    "max_AMC_nut_reduction", "max_SRSA_nut_reduction"];
inf_p = (; zip(Symbol.(param_names), param_vals)...)

###########################

let
    plotID = "SEG12"
    data, sol = get_plottingdata(sim;
        plotID,
        inf_p,
        nspecies = 25,
        startyear = 2009,
        endyear = 2021,
        seed = rand(1:1000))

    selected_patch = 1
    simbiomass = TimeArray((;
            biomass = vec(sum(ustrip.(sol.biomass[:, selected_patch, :]); dims = 2)),
            numeric_date = sol.numeric_date,
            date = sol.date), timestamp = :date)

    mbiomass = data.measured_biomass

    simbiomass_sub = simbiomass[timestamp(data.measured_biomass)]
    simbiomass_vals = values(simbiomass_sub.biomass)

    if any(isnan.(simbiomass_vals))
        @warn "Biomass NaN"
        return -Inf
    end

    ### calculate the likelihood
    biomass_d = Laplace.(simbiomass_vals, inf_p.sigma_biomass)
    ll_measuredbiomass = logpdf.(biomass_d, values(data.measured_biomass.biomass))

    fig = Figure()
    Axis(fig[1, 1]; title = plotID)
    lines!(values(simbiomass.numeric_date), values(simbiomass.biomass);
        color = :orange)
    scatter!(values(mbiomass.numeric_date), values(mbiomass.biomass);
        color = :black, markersize = 8)
    text!(values(mbiomass.numeric_date), values(mbiomass.biomass);
        text = string.(round.(ll_measuredbiomass; digits = 1)))
    fig
end
