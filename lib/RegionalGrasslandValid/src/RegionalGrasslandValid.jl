module RegionalGrasslandValid

import CSV
import Dates
using DataFrames, DataFramesMeta
using Unitful
using Distributions
using LinearAlgebra

include("input_data.jl")
include("validation_data.jl")
include("ll_calculation.jl")
include("optim_calculation.jl")
include("traits/Traits.jl")

export loglikelihood_model, get_plottingdata, get_validation_data
export ll_VIPS_t, ll_VIPS
export optimize_biomass

function __init__()
    @info "Loading data of RegionalGrasslandValid"
    datapath = "../lib/RegionalGrasslandData"
    load_data(datapath)
    Traits.load_data(datapath)
end

function load_data(datapath)
    ########### validation data
    soilmoisture = CSV.read("$datapath/validation/soilmoisture.csv",
        DataFrame)

    evaporation = CSV.read("$datapath/validation/evaporation.csv",
        DataFrame)

    satbiomass = CSV.read("$datapath/validation/sat_biomass.csv",
        DataFrame)

    measuredveg = CSV.read("$datapath/validation/measured_veg.csv",
        DataFrame)

    valid = (;
        soilmoisture,
        evaporation,
        satbiomass,
        measuredveg)

    ########### input data
    initbiomass = CSV.read("$datapath/input/init_biomass.csv",
        DataFrame)

    clim = CSV.read("$datapath/input/temperature_precipitation.csv",
        DataFrame)

    pet = CSV.read("$datapath/input/PET.csv",
        DataFrame)

    par = CSV.read("$datapath/input/par.csv",
        DataFrame)

    ### mean index from 2011, 2014, 20117, 2021
    nut = CSV.read("$datapath/input/soilnutrients.csv",
        DataFrame)

    soil = CSV.read("$datapath/input/soilwater.csv",
        DataFrame)

    mow = CSV.read("$datapath/input/mowing.csv",
        DataFrame)

    graz = CSV.read("$datapath/input/grazing.csv",
        DataFrame)

    input = (;
        initbiomass,
        clim,
        pet,
        par,
        nut,
        soil,
        mow,
        graz)

    global data = (;
        input,
        valid)

    return nothing
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
    sol = sim.solve_prob(; input_obj)

    ########################## Measured data
    data = get_validation_data(; plotID)

    return data, sol
end

end
