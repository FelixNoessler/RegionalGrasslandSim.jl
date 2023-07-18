module RegionalGrasslandValid

include("traits/Traits.jl")

import CSV
import Dates
using DataFrames, DataFramesMeta
using Unitful
using Distributions
using LinearAlgebra

export loglikelihood_model, get_plottingdata, get_validation_data
export ll_VIPS_t, ll_VIPS


struct DataTables
    soilwater
    satbiomass
    measuredveg
    clim
    pet
    par
    nut
    soil
    mow
    graz
end

function __init__()
    @info "Loading data of RegionalGrasslandValid"
    datapath = "../RegionalGrasslandData"
    load_data(datapath)
    Traits.load_data(datapath)
end

VIP_plots = ["$(explo)0$i" for i in 1:9 for explo in ["HEG", "SEG", "AEG"]];

function load_data(datapath)

    ########### validation data
    soilwater_df = CSV.read(
        "$datapath/validation/soilwater.csv",
        DataFrame)

    satbiomass_df = CSV.read(
        "$datapath/validation/sat_biomass.csv",
        DataFrame
    )

    measuredveg_df = CSV.read(
        "$datapath/validation/measured_veg.csv",
        DataFrame
    )

    ########### input data
    clim_df = CSV.read(
        "$datapath/input/climate.csv",
        DataFrame)

    pet_df = CSV.read(
        "$datapath/input/PET.csv",
        DataFrame
    )

    par_df = CSV.read(
        "$datapath/input/par.csv",
        DataFrame
    )

    ### mean index from 2011, 2014, 20117, 2021
    nut_df = CSV.read(
        "$datapath/input/soilnutrients.csv",
        DataFrame
    )

    soil_df = CSV.read(
        "$datapath/input/soilwater.csv",
        DataFrame
    )

    mow_df = CSV.read(
        "$datapath/input/mowing.csv",
        DataFrame
    )

    graz_df = CSV.read(
        "$datapath/input/grazing.csv",
        DataFrame
    )

    global dat = DataTables(
        soilwater_df,
        satbiomass_df,
        measuredveg_df,
        clim_df,
        pet_df,
        par_df,
        nut_df,
        soil_df,
        mow_df,
        graz_df
    )

    return nothing
end


function ll_VIPS_t(sim; inf_p)
    ll = Threads.Atomic{Float64}(0.0)
    Threads.@threads for plotID in VIP_plots
        ll_plot = loglikelihood_model(sim;
            plotID,
            inf_p,
            startyear=2012,
            endyear=2021)

        Threads.atomic_add!(ll, ll_plot)
    end

    ### free RAM space
    GC.gc()
    ccall(:malloc_trim, Cvoid, (Cint,), 0)

    return ll[]
end

function ll_VIPS(sim; inf_p)
    ll = 0.0
    for plotID in VIP_plots
        ll_plot = loglikelihood_model(sim;
            plotID,
            inf_p,
            startyear=2012,
            endyear=2021)

        ll += ll_plot
    end

    return ll
end


function yearly_temp_cumsum(d)
    adj = 365
    nyears = length(d) ÷ adj + 1

    final_cumsum = Array{Float64}(undef, length(d))

    d = ustrip.(d)
    d[d .< 0.0] .= 0.0

    for y in 1:nyears
        sliced_d = d[1+adj*(y-1):min(adj*y, length(d))]
        final_cumsum[1+adj*(y-1):min(adj*y, length(d))] .= cumsum(sliced_d)
    end

    return final_cumsum
end

function validation_input(;
        plotID,
        inf_p,
        nyears)

    nspecies = 25
    npatches = 1
    water_reduction = true
    nutrient_reduction = true

    clim_sub = @subset dat.clim :plotID .== plotID
    pet_sub = @subset dat.pet first.(:explo) .== first(plotID)
    par_sub = @subset dat.par first.(:explo) .== first(plotID)
    nut_sub = @subset dat.nut :plotID .== plotID
    soil_sub = @subset dat.soil :plotID .== plotID
    mow_sub = @subset dat.mow :plotID .== plotID
    graz_sub = @subset dat.graz :plotID .== plotID

    ### ----------------- abiotic
    temperature = clim_sub.Ta_200 .* u"°C"
    temperature_sum = yearly_temp_cumsum(temperature)
    precipitation = clim_sub.precipitation .* u"mm / d"
    PAR = par_sub.PAR .* u"MJ / (d * m^2)"
    PET = pet_sub.PET .* u"mm / d"
    nutrient_index = nut_sub.n_index[1]
    WHC = soil_sub.WHC[1] * u"mm"
    PWP = soil_sub.PWP[1] * u"mm"
    Clay = soil_sub.Clay[1]
    Silt = soil_sub.Silt[1]
    Sand = soil_sub.Sand[1]
    organic = soil_sub.organic[1]
    bulk = soil_sub.bulk[1]
    root_depth = soil_sub.root_depth[1]

    ### ----------------- traits
    traits = Traits.random_traits(nspecies;);
    sort!(traits, :SLA)
    relative_traits = Traits.relative_traits(; trait_data=traits)

    ### ----------------- mowing
    mowing_day = [mow_sub[year, "MowingDay$i"] for year in 1:nyears, i in 1:5]

    mowing_height = [mow_sub[year, "CutHeight_cm$i"] for year in 1:nyears, i in 1:5]
    height_missing = ismissing.(mowing_height) .&& .! ismissing.(mowing_day)


    mowing_height = convert(Matrix{Union{Missing, Int64}}, mowing_height)
    mowing_height[height_missing] .= 7
    height_should_be_missing = .! ismissing.(mowing_height) .&& ismissing.(mowing_day)
    mowing_height[height_should_be_missing] .= missing

    mowing_heights = [mowing_height[y, :][
        .! ismissing.(mowing_height[y, :])
    ] for y in 1:nyears]
    mowing_heights = convert(Vector{Vector{Int64}}, mowing_heights)

    mowing_days = [mowing_day[y, :][.! ismissing.(mowing_day[y, :])] for y in 1:nyears]
    mowing_days = convert(Vector{Vector{Int64}}, mowing_days)

    ### ----------------- grazing
    grazing_start = [
        graz_sub[year, "StartGrazingPeriod$i"] for year in 1:nyears, i in 1:4
    ]
    days_grazing = [
        graz_sub[year, "DayGrazing$i"] for year in 1:nyears, i in 1:4
    ]

    ### set the start day to missing, if grazing period == 0
    grazing_start = convert(Matrix{Union{Missing, Int64}}, grazing_start)
    grazing_start[iszero.(days_grazing)] .= missing

    ### calculate the end day of the grazing period
    grazing_end = grazing_start .+ days_grazing

    ### get the grazing intensity and transform 0 and NaNs to missing
    grazing_intensity = [
        graz_sub[year, "GrazingIntensity$i"] for year in 1:nyears, i in 1:4
    ]
    grazing_intensity = convert(Matrix{Union{Missing, Float64}}, grazing_intensity)
    intensity_filter = isnan.(grazing_intensity) .|| iszero.(grazing_intensity)
    grazing_intensity[intensity_filter] .= missing

    ### convert matrices with missing to vectors of vectors
    grazing_start = [
        grazing_start[y, :][.! ismissing.(grazing_start[y, :])] for y in 1:nyears
    ]
    grazing_end = [
        grazing_end[y, :][.! ismissing.(grazing_end[y, :])] for y in 1:nyears
    ]
    grazing_intensity = [
        grazing_intensity[y, :][.! ismissing.(grazing_intensity[y, :])] for y in 1:nyears
    ]

    ### derive the final grazing density vector
    grazing_density = fill(0.0u"ha ^ -1", 365*nyears)
    for y in 1:nyears
        for i in eachindex(grazing_start[y])
            start_day = 365 * (y-1) + grazing_start[y][i]
            end_day = 365 * (y-1) + grazing_end[y][i]
            livestock_density = grazing_intensity[y][i]

            if isa(start_day, Number) && isa(end_day, Number) && isa(livestock_density, Number)
                end_day = min(365 * nyears, end_day)

                grazing_density[start_day : end_day] .=
                    livestock_density .* u"ha ^ -1"
            end
        end
    end

    return (
        inf_p,
        env_data = (;
            PAR,
            precipitation,
            temperature,
            temperature_sum,
            PET,
        ),
        site=(;
            nutrient_index,
            WHC,
            PWP,
            Clay,
            Silt,
            Sand,
            organic,
            bulk,
            root_depth
        ),
        traits,
        relative_traits,
        mowing_days,
        mowing_heights,
        grazing=grazing_density,
        nspecies,
        npatches,
        water_reduction,
        nutrient_reduction,
        nyears
    )
end



function get_validation_data(; plotID)
    soilwater_sub = @subset dat.soilwater :plotID .== plotID
    satbiomass_sub = @subset dat.satbiomass :plotID .== plotID
    measuredveg_sub = @subset dat.measuredveg :plotID .== plotID

    return (;
        soil_moisture=soilwater_sub.soil_moisture,
        evaporation=soilwater_sub.evaporation,
        biomass=satbiomass_sub.biomass_kg_ha,
        biomass_t=satbiomass_sub.t,
        measured_biomass=measuredveg_sub.biomass_kg_ha,
        measured_biomass1=measuredveg_sub.biomass1_kg_ha,
        measured_biomass_t=measuredveg_sub.t,
    )
end


function get_plottingdata(sim::Module;
    inf_p,
    plotID,
    startyear,
    endyear)

    ########################## Run model
    nyears = length(startyear:endyear)
    input_obj = validation_input(;
        plotID, nyears, inf_p)
    sol = sim.solve_prob(; input_obj);

    ########################## Measured data
    data = get_validation_data(; plotID)

    return data, sol
end


function loglikelihood_model(sim::Module;
        inf_p,
        plotID,
        startyear=2012,
        endyear=2021,
        inityears=5)

    inf_p = (; inf_p ...)
    data, sol = get_plottingdata(sim;
        inf_p,
        plotID,
        startyear,
        endyear)

    ########################## Calculate likelihood
    ###### biomass
    biomass_sum = vec(sum(ustrip.(sol.biomass); dims=3))[2:end][1+365*inityears:end]

    if any(isnan.(biomass_sum))
        @error "Biomass sum isnan"
    end

    # if iszero(biomass_sum[end])
    #     return -Inf
    # end

    # f = 0 .< data.measured_biomass_t .< 1825
    # sim_biomass = biomass_sum[data.measured_biomass_t[f]]
    # biomass_d = MvNormal(sim_biomass, inf_p.sigma_biomass * I)
    # ll_biomass = logpdf(biomass_d, data.measured_biomass[f])

    sim_biomass = biomass_sum[data.biomass_t]
    biomass_d = MvNormal(sim_biomass, inf_p.sigma_biomass * I)
    ll_biomass = logpdf(biomass_d, data.biomass)

    ###### evaporation
    ll_evaporation = 0.0
    evapo_index = .! ismissing.(data.evaporation)
    if ! all(iszero.(evapo_index))
        sim_evaporation = ustrip.(sol.evaporation[2:end])[1+365*inityears:end][evapo_index]
        data_evaporation = float.(data.evaporation[evapo_index])
        evaporation_d = MvNormal(sim_evaporation, inf_p.sigma_evaporation * I)
        ll_evaporation = logpdf(evaporation_d, data_evaporation)
    end

    # ###### soil moisture
    sim_soilmoisture = ustrip.(sol.water)[2:end][1+365*inityears:end] ./ sol.p.site.root_depth .* inf_p.moisture_conv
    soilmoisture_d = MvNormal(sim_soilmoisture, inf_p.sigma_soilmoisture * I)
    ll_soilmoisture = logpdf(soilmoisture_d, data.soil_moisture)

    ###### total log likelihood
    ll = ll_biomass #+ ll_evaporation + ll_soilmoisture

    return ll
end

end
