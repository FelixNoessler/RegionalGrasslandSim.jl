# ] dev ../
# ] dev ../lib/RegionalGrasslandVis/
# ] dev ../lib/RegionalGrasslandValid/
# ] dev ../lib/RegionalGrasslandScenario/

using Unitful
using Statistics

import RegionalGrasslandSim as sim
import RegionalGrasslandScenario as scen
import RegionalGrasslandValid as valid
import RegionalGrasslandVis as vis
using RegionalGrasslandValid

#############################
# using JuliaFormatter, RegionalGrasslandSim
# format(joinpath(dirname(pathof(RegionalGrasslandSim)), ".."))

#############################
param_names = [
    "moistureconv_alpha", "moistureconv_beta",
    "senescence_intercept", "senescence_rate",
    "below_competition_strength", "trampling_factor", "grazing_half_factor",
    "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
    "max_AMC_nut_reduction", "max_SRSA_nut_reduction",
    "b_biomass",
    "b_SLA", "b_LNCM", "b_AMC", "b_height", "b_SRSA_above",
    "b_soilmoisture"];
vals = [
    28.428282890792907, 72.70570950774605, 4.999073781961289, 4.153059795872836,
    0.3931281548093227, 267.95272883632566, 1874.387202377206, 10.23688734649984,
    0.44395619127597463, 0.2388901001112616, 0.36320671631375023, 0.5711510823547715,
    1628.6040232644204, 2145.2018690252116, 1021.0380327120577, 8.80738136284001,
    260.21386563130073, 2129.699173441751, 87.03525031931032]
inf_p = (; zip(Symbol.(param_names), vals)...)
input_obj = valid.validation_input(;
    plotID = "HEG01", nspecies = 25,
    startyear = 2009, endyear = 2021,
    inf_p, npatches = 1);
@time sol = sim.solve_prob(; input_obj);

# @time data, sol = get_plottingdata(sim;
#     plotID,
#     inf_p,
#     nspecies = 25,
#     startyear = 2009,
#     endyear = 2021,
#     seed = rand(1:100));

training_plots = ["$(explo)$(lpad(i, 2, "0"))" for i in 1:9
                  for explo in ["HEG", "SEG", "AEG"]]

sum([loglikelihood_model(sim;
    inf_p,
    plotID = p,
    nspecies = 30) for p in training_plots])

############################# Dashboard
using GLMakie
GLMakie.activate!()
Makie.inline!(false)
vis.dashboard(; sim, valid, scen, inf_p_start = inf_p)

############################# Validation
using CairoMakie
Makie.inline!(true)

let
    plotID = "HEG01"
    trait = :height
    patch = 1
    fig = Figure()
    Axis(fig[1, 1])

    for i in [1]
        input_obj = RegionalGrasslandValid.validation_input(;
            plotID, nspecies = 50,
            startyear = 2009, endyear = 2021,
            inf_p)
        sol = sim.solve_prob(; input_obj)
        scatter!(ustrip.(sol.p.species[trait]),
            ustrip.(sol.biomass[end, patch, :]);
            color = :red)
    end
    fig
end

let
    plotID = "SEG01"
    inf_p = (; zip(Symbol.(param_names), param_vals)...)
    vis.biomass_validation2(sim, valid; plotID, inf_p)
end

let
    # selected_plots = ["$(explo)$(lpad(i, 2, "0"))" for i in 1:50
    #                   for explo in ["HEG", "SEG", "AEG"]]

    selected_plots = ["SEG07"]
    inf_p = (; zip(Symbol.(param_names), param_vals)...)

    for plotID in selected_plots
        @info plotID

        # loglikelihood_model(sim;
        #     inf_p,
        #     plotID,
        #     nspecies=25)
        # data, sol = get_plottingdata(sim;
        #     plotID,
        #     inf_p,
        #     startyear = 2012,
        #     endyear = 2021)

        # f = vis.biomass_validation(data, sol; plotID)
        # save("../tmp/0_$plotID.png", f)

        # f = vis.soilmoisture_validation(data, sol;
        #     converted = true)
        # save("../tmp/1_water_$plotID.png", f)

        # f = vis.evaporation_validation(data, sol)
        # save("../tmp/2_evapo_$plotID.png", f)
    end
end

param_vals = [
    35914.220404412874,
    70.43347155395173,
    75.00103265678999,
    0.2114621486366593,
    0.4538043668775889,
    32.78528981780835,
    3.865860365279632,
    174.25334577381417,
    463.81007069310135,
    15.633385764983174,
    0.3499484166680329,
    0.41143613204926305,
    0.5788642998287264,
    0.5981211370597104,
]
param_names = [
    "sigma_biomass", "sigma_evaporation", "sigma_soilmoisture",
    "moisture_conv", "senescence_intercept", "senescence_rate",
    "below_competition_strength", "trampling_factor", "grazing_half_factor",
    "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
    "max_AMC_nut_reduction", "max_SRSA_nut_reduction"]
inf_p = (; zip(Symbol.(param_names), param_vals)...)

@time data, sol = get_plottingdata(sim;
    plotID = "HEG01",
    inf_p,
    startyear = 2012,
    endyear = 2021);

input_obj = valid.validation_input(;
    plotID = "HEG01",
    inf_p,
    nyears = length(2012:2021))

@time valid.Traits.calc_CWM(;
    biomass = sol.biomass,
    trait_data = sol.p.species.SLA);
@time valid.Traits.calc_CWV(;
    biomass = sol.biomass,
    trait_data = sol.p.species.SLA);

############################# Run one simulation
# using CairoMakie
# CairoMakie.activate!()
# let
#     nspecies = 5;
#     npatches = 1;
#     nyears = 5;
#     explo = "HAI";
#     input_obj = prepare_data(; nspecies, npatches, nyears, explo)

#     ########################## run model
#     @time sol = sim.solve_prob(;
#         input_obj);
#     nothing
#     # vis.mowing(sim;
#     #     nspecies,
#     #     λ=sol.p.species.λ)

#     # vis.trampling(sim;)

#     # vis.trampling_combined(sim;
#     #     LA=input_obj.traits.LA,
#     #     height=input_obj.traits.height,
#     #     nspecies)
#     # vis.seasonal_component_senescence(sim;
#     #     STs=input_obj.env_data.temperature_sum)
#     # vis.seasonal_effect(sim;
#     #     STs=input_obj.env_data.temperature_sum)
#     # vis.temperatur_reducer(sim;
#     #     Ts=input_obj.env_data.temperature)
#     # vis.radiation_reducer(sim;
#     #     PARs=input_obj.env_data.PAR)

#     # vis.climate_series(; sol)

#     # vis.sla_water_response(sim;
#     #     nspecies,
#     #     SLA=input_obj.traits.SLA)

#     # vis.potential_growth(sim;
#     #     nspecies,
#     #     SLA=input_obj.traits.SLA
#     # )

#     # vis.myco_response(sim;
#     #     nspecies,
#     #     mycorrhizal_colon=input_obj.traits.AMC
#     # )

#     # vis.srsa_response(sim;
#     #     nspecies,
#     #     SRSA_above=input_obj.traits.SRSA_above)

#     # sol.biomass

#     # sol.t

#     # sol.p.site
# end
