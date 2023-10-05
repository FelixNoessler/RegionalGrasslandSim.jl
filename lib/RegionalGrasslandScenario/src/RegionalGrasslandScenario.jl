module RegionalGrasslandScenario

using Unitful
using DataFrames
using TimeSeries

include("climate_analysis/AirTemperature.jl")
include("climate_analysis/PAR.jl")
include("climate_analysis/Precipitation.jl")
include("climate_analysis/PET.jl")
include("landuse/Landuse.jl")

export scenario_input

function __init__()
    @info "Loading data of RegionalGrasslandScenario"
    datapath = "../lib/RegionalGrasslandData"

    AirTemperature.load_data(datapath)
    PAR.load_data(datapath)
    Precipitation.load_data(datapath)
    PET.load_data(datapath)

    return nothing
end

function scenario_input(;
    nyears,
    nspecies,
    explo,
    mowing_doys,
    grazing_start,
    grazing_end,
    grazing_intensity,
    nutrient_index,
    WHC,
    PWP,
    npatches = 1,
    constant_seed = false,
    initbiomass = 500u"kg/ha",
    senescence_included = true,
    potgrowth_included = true,
    mowing_included = true,
    grazing_included = true,
    below_included = true,
    height_included = true,
    water_red = true,
    nutrient_red = true,
    temperature_red = true,
    season_red = true,
    radiation_red = true)
    temp_data = AirTemperature.predict_temperature(;
        nyears,
        explo)
    par_data = PAR.predict_par(;
        explo,
        nyears)
    precipitation_data = Precipitation.predict_precipitation(;
        nyears,
        explo)
    evapo_data = PET.predict_pet(;
        explo,
        nyears)
    grazing_data = Landuse.grazing_input(;
        grazing_start,
        grazing_end,
        grazing_intensity,
        nyears)

    mowing_data = Landuse.mowing_input(;
        mowing_doys,
        nyears,
        cutting_height = 0.07u"m")

    d = Dates.Date(0):Dates.lastdayofyear(Dates.Date(nyears - 1))
    f = Dates.dayofyear.(d) .<= 365

    daily_data = (;
        date = d[f],
        temperature = temp_data,
        temperature_sum = AirTemperature.yearly_temp_cumsum(temp_data),
        precipitation = precipitation_data,
        PET = evapo_data,
        PAR = par_data,
        mowing = mowing_data,
        grazing = grazing_data)

    #### -------------- whether parts of the simulation are included
    included = (;
        senescence_included,
        potgrowth_included,
        mowing_included,
        grazing_included,
        below_included,
        height_included,
        water_red,
        nutrient_red,
        temperature_red,
        season_red,
        radiation_red)

    return (
        nspecies,
        npatches,
        patch_xdim = Int(sqrt(npatches)),
        patch_ydim = Int(sqrt(npatches)),
        constant_seed,
        included,
        startyear = 0,
        endyear = nyears,
        site = (;
            nutrient_index,
            WHC = WHC * u"mm",
            PWP = PWP * u"mm",
            initbiomass),
        daily_data)
end

end # module RegionalGrasslandScenario
