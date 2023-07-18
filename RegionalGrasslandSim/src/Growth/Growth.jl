module Growth

using Unitful

include("growth reducers.jl")
include("defoliation.jl")
include("senescence.jl")

"""
    growth(; nspecies, fun_response, SLA, biomass, PAR, WR, T, ST)

TBW
"""
function growth(;
        nspecies,
        fun_response,
        trait_similarity,
        below_competition_strength,
        SLA, CH, biomass,
        PAR,
        WR,
        PWP, WHC,
        nutrients,
        PET,
        T, ST,
        water_red,
        nutrient_red,
        max_SRSA_water_reduction,
        max_SLA_water_reduction,
        max_AMC_nut_reduction,
        max_SRSA_nut_reduction)

    #### potential growth
    pot_growth, lai_tot = potential_growth(; nspecies, SLA, biomass, PAR)

    ### influence of the height of plants
    cwm_CH = community_weighted_mean_height(; biomass, CH)
    CHinfluence = (CH .- (CH .- cwm_CH) ./ 1.5 ) ./ cwm_CH

    #### below ground competition --> trait similarity and abundance
    below = below_ground_competition(;
        biomass, trait_similarity, nspecies,
        below_competition_strength)

    #### growth reducer
    Rred = radiation_reduction(; PAR)
    Tred = temperature_reduction(; T)
    Waterred = water_reduction(;
        fun_response, WR, water_red, PET,
        max_SRSA_water_reduction,
        max_SLA_water_reduction,
        PWP, WHC)
    Nutred = nutrient_reduction(;
        fun_response, nutrient_red,
        nutrients,
        max_AMC_nut_reduction,
        max_SRSA_nut_reduction)
    Seasonalred = seasonal_reduction(; ST)

    reduction = Rred .* Tred .* Waterred .* Nutred .* Seasonalred .* CHinfluence .* below

    #### final growth
    return pot_growth .* reduction, lai_tot
end


function below_ground_competition(;
    biomass, trait_similarity, nspecies,
    below_competition_strength)

    reduction_coefficient = Array{Float64}(undef, nspecies)
    biomass = ustrip.(biomass)

    for i in 1:nspecies
        x = sum(trait_similarity[i, :] .* biomass)
        reduction_coefficient[i] = exp(-below_competition_strength * x)
    end

    return reduction_coefficient
end

@doc raw"""
    potential_growth(; SLA, nspecies, biomass, PAR)

Calculates the potential growth of all plant species
in a specific patch.

This function is called each time step (day) for each patch.
The NamedTuple `p` contains all the species specific trait values.
The vector `biomass` contains the biomass of the species in
the specific patch. The `PAR` value is the photosynthetically
active radiation of the day.

First, the leaf area indices of all species are calculated
(see [`calculate_LAI`](@ref)). Then, the total leaf area is
computed. An inverse exponential function is used to calculate
the total primary production:

```math
\text{totalgrowth} = 10 \cdot PAR \cdot RUE_{max} \cdot (1 -  \text{exp}(-\alpha \cdot \text{LAItot}))
```

This primary production is then multiplied with the share of the
leaf area index of the individual species

![Influence of the specific leaf area on the potential growth](../../img/sla_potential_growth.svg)
"""
function potential_growth(; SLA, nspecies, biomass, PAR)
    lais = calculate_LAI(; SLA, biomass)
    LAItot = sum(lais)

    if LAItot == 0
        return fill(0.0, nspecies)u"kg / ha / d", LAItot
    end

    RUE_max = 3u"g / MJ" # Maximum radiation use efficiency g DM MJ-1
    α = 0.6   # Extinction coefficient, unitless

    total_growth = PAR * RUE_max * (1 - exp(-α * LAItot))
    species_growth = total_growth .* lais ./ LAItot

    return uconvert.(u"kg / ha / d", species_growth), LAItot
end


"""
    community_weighted_mean_height(; biomass, CH)

TBW
"""
function community_weighted_mean_height(; biomass, CH)
    return sum(biomass .* CH) / sum(biomass)
end


@doc raw"""
    calculate_LAI(; SLA, biomass)

Calculate the leaf area index of all species of one habitat patch.
"""
function calculate_LAI(; SLA, biomass)
    LAM = 0.62 # Proportion of laminae in green biomass
    return uconvert.(NoUnits, SLA .* biomass .* LAM)
end

end # of module
