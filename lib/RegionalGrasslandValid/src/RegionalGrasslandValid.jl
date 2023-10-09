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

export loglikelihood_model, get_plottingdata, get_validation_data
export ll_VIPS_t, ll_VIPS
export optimize_biomass

function __init__()
    @info "Loading data of RegionalGrasslandValid"
    datapath = "../lib/RegionalGrasslandData"
    load_data(datapath)
    model_parameters()
end

struct ModelParam
    names::Vector{String}
    best::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
end

function model_parameters()
    names = [
        "moistureconv_alpha", "moistureconv_beta",
        "senescence_intercept", "senescence_rate",
        "belowground_density_effect", "belowtrait_similarity_exponent",
        "height_strength", "leafnitrogen_graz_exp",
        "trampling_factor", "grazing_half_factor",
        "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
        "max_AMC_nut_reduction", "max_SRSA_nut_reduction",
        "b_biomass",
        "b_SLA", "b_LNCM", "b_AMC", "b_height", "b_SRSA_above",
        "b_soilmoisture"]
    best = [
        15.974238060436035,
        283.8284074247446,
        6.58833885841744,
        2.680997367833937,
        1.0346464525351673,
        19.71635483939154,
        0.5,
        3,
        251.43022955138898,
        1645.3988036165588,
        3.523324632368246,
        0.7561526618457965,
        0.8555098929939224,
        0.00013975735455720893,
        0.6368009957917703,
        1238.3915632978947,
        4069.5097732255103,
        32.21560656121319,
        316.9027469361824,
        3638.0439110435286,
        113.86586008943979,
        405.40235861163285,
    ]

    #         mc_α mc_β s_i s_r below bsim height ni_graz tram graz mow SRSA SLA  AMC  SRSA_n
    lb_prep = [0,    0,  0,  0, -10,   0,  0,     0,      100,    0,  0, 0.0, 0.0, 0.0, 0.0]
    ub_prep = [80, 300, 10, 10, 2.5,   30, 1,     15,     300, 2000, 50, 1.0, 1.0, 1.0, 1.0]
    nscale_params = sum(startswith.(names, "b_"))
    lb_b = zeros(nscale_params)
    ub_b = fill(5e3, nscale_params)
    lb = vcat(lb_prep, lb_b)
    ub = vcat(ub_prep, ub_b)

    return ModelParam(names, best, lb, ub)
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
    input_objs,
    inf_p,
    plotID,
    startyear)

    ########################## Run model
    sol = sim.solve_prob(; input_obj = input_objs[plotID], inf_p)

    ########################## Measured data
    ## I shouldn't call this function each time...
    data = get_validation_data(; plotID, startyear)

    return data, sol
end

end
