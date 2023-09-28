module RegionalGrasslandValid

import CSV
import Dates
using DataFrames, DataFramesMeta
using Unitful
using Distributions
using LinearAlgebra
using TimeSeries

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

    measuredbiomass = CSV.read("$datapath/validation/measured_biomass.csv",
        DataFrame)

    mtraits = CSV.read("$datapath/validation/cwm_traits.csv",
        DataFrame)

    traits = (vals = [mtraits.SLA mtraits.LNCM mtraits.AMC mtraits.SRSA_above mtraits.height],
        dim = [:SLA, :LNCM, :AMC, :SRSA_above, :height],
        t = mtraits.date,
        num_t = mtraits.numeric_date,
        plotID = mtraits.plotID)

    valid = (;
        soilmoisture,
        evaporation,
        traits,
        measuredbiomass)

    ########### input data
    initbiomass = CSV.read("$datapath/input/init_biomass.csv",
        DataFrame)

    ## time dependent 2009-2022
    clim = CSV.read("$datapath/input/temperature_precipitation.csv",
        DataFrame)

    ## time dependent 2006-2022
    pet = CSV.read("$datapath/input/PET.csv",
        DataFrame)

    ## time dependent 2006-2022
    par = CSV.read("$datapath/input/par.csv",
        DataFrame)

    ### mean index from 2011, 2014, 20117, 2021
    nut = CSV.read("$datapath/input/soilnutrients.csv",
        DataFrame)

    ## constant WHC & PWP
    soil = CSV.read("$datapath/input/soilwater.csv",
        DataFrame)

    ## time dependent 2006 - 2021
    mow = CSV.read("$datapath/input/mowing.csv",
        DataFrame)

    ## time dependent 2006 - 2021
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
    nspecies,
    startyear,
    endyear,
    seed)

    ########################## Run model
    input_obj = validation_input(;
        plotID, nspecies, startyear, endyear, inf_p, seed)
    sol = sim.solve_prob(; input_obj)

    ########################## Measured data
    data = get_validation_data(; plotID, startyear)

    return data, sol
end

end
