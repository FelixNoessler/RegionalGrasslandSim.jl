module RegionalGrasslandScenario

using Unitful
using DataFrames

include("climate_analysis/AirTemperature.jl")
include("climate_analysis/PAR.jl")
include("climate_analysis/Precipitation.jl")
include("climate_analysis/PET.jl")
include("landuse/Landuse.jl")
include("traits/Traits.jl")

export scenario_input

function __init__()
    @info "Loading data of RegionalGrasslandScenario"
    datapath = "../RegionalGrasslandData"

    AirTemperature.load_data(datapath)
    PAR.load_data(datapath)
    Precipitation.load_data(datapath)
    PET.load_data(datapath)
    Traits.load_data(datapath)

    return nothing
end

function scenario_input(;
        explo, nyears, nspecies,
        inf_p,
        npatches=1,
        mowing_days=[Vector{Int64}() for _ in 1:nyears],
        mowing_heights=[Vector{Int64}() for _ in 1:nyears],
        nutrient_index,
        WHC,
        PWP,
        grazing_start=[],
        grazing_end=[],
        grazing_intensity=[],
        water_reduction,
        nutrient_reduction)

    temp_data = AirTemperature.predict_temperature(;
        nyears,
        explo);
    par_data = PAR.predict_par(;
        explo,
        nyears);
    precipitation_data = Precipitation.predict_precipitation(;
        nyears,
        explo)
    evapo_data = PET.predict_pet(;
        explo,
        nyears);
    trait_data = Traits.random_traits(nspecies;);
    relative_trait_data = Traits.relative_traits(;trait_data)
    grazing_data = Landuse.grazing_input(;
        grazing_start,
        grazing_end,
        grazing_intensity,
        nyears
    )

    return (
        env_data = (;
            PAR=par_data,
            precipitation=precipitation_data,
            temperature=temp_data,
            temperature_sum=AirTemperature.yearly_temp_cumsum(temp_data),
            PET=evapo_data,
        ),
        inf_p,
        site=(;
            nutrient_index,
            WHC=WHC * u"mm",
            PWP=PWP * u"mm"
        ),
        traits=trait_data,
        relative_traits=relative_trait_data,
        mowing_days,
        mowing_heights,
        grazing=grazing_data,
        nutrient_index,
        nspecies,
        npatches,
        water_reduction,
        nutrient_reduction,
        nyears
    )
end

end # module RegionalGrasslandScenario
