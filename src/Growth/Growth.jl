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
    CHinfluence = (CH .- (CH .- cwm_CH) ./ 1.5) ./ cwm_CH

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

@doc raw"""
    similarity_matrix(; scaled_traits)

Computes the trait similarity of all plant species.

The trait similarity between plant species $i$ and
plant species $u$ for $T$ traits is calculated as follows:
```math
\text{trait_similarity}_{i,u} =
    1-\frac{\sum_{t=1}^{t=T}
        |\text{scaled_trait}_{t,i} - \text{scaled_trait}_{t,u}|}{T}
```

To give each functional trait an equal influence,
the trait values have been scaled by the 5 % ($Q_{0.05, t}$)
and 95 % quantile ($Q_{0.95, t}$) of trait values of 100 plant species:
```math
\text{scaled_trait}_{t,i} =
    \frac{\text{trait}_{t,i} - Q_{0.05, t}}
    {Q_{0.95, t} - Q_{0.05, t}}
```

If the rescaled trait values were below zero or above one, the values were
set to zero or one respectively.
"""
function similarity_matrix(; scaled_traits)
    nspecies, ntraits = size(scaled_traits)
    similarity_mat = Array{Float64}(undef, nspecies, nspecies)

    for i in 1:nspecies
        diff = zeros(nspecies)

        for n in 1:ntraits
            species_trait = scaled_traits[i, n]
            diff .+= abs.(species_trait .- scaled_traits[:, n])
        end

        similarity_mat[i, :] .= 1 .- (diff ./ ntraits)
    end

    return similarity_mat
end

@doc raw"""
    below_ground_competition(;
        biomass, trait_similarity, nspecies,
        below_competition_strength)

Models the below-ground competiton between plant.

Plant growth is reduced if a large biomass of plant species with similar
functional traits is already present. The `below_competition` factor has
a value between 0 and 1. For plant species $i$ with $N$ plant species present
it is defined as follows:

```math
\text{below_competition}_i =
    exp\left(-\frac{\text{below_competition_strength}}{1000} \cdot
        \left[\sum_{u=1}^{u=N} \text{trait_similarity}_{i,u} \cdot \text{biomass}_u\right]
    \right)
```

The `below_competition_strength` can therefore be seen as a parameter that
controls the density dependence.

The `trait_similarity` is computed by the function [`similarity_matrix`](@ref)
and includes the traits specific leaf area (`SLA`),
arbuscular mycorrhizal colonisation rate (`AMC`),
and the root surface area devided by the above ground biomass (`SRSA_above`).
"""
function below_ground_competition(;
    biomass, trait_similarity, nspecies,
    below_competition_strength)
    reduction_coefficient = Array{Float64}(undef, nspecies)

    for i in 1:nspecies
        x = sum(trait_similarity[i, :] .* biomass)
        reduction_coefficient[i] = exp(-below_competition_strength / 1000 * x)
    end

    return reduction_coefficient
end

@doc raw"""
    potential_growth(; SLA, nspecies, biomass, PAR)

Calculates the potential growth of all plant species
in a specific patch.

This function is called each time step (day) for each patch.
The `PAR` value is the photosynthetically
active radiation of the day.

First, the leaf area indices of all species are calculated
(see [`calculate_LAI`](@ref)). Then, the total leaf area is
computed. An inverse exponential function is used to calculate
the total primary production:

```math
\text{totalgrowth} =
    PAR \cdot RUE_{max} \cdot (1 -  \text{exp}(-\alpha \cdot \text{LAItot}))
```

This primary production is then multiplied with the share of the
leaf area index of the individual species

![Influence of the specific leaf area on the potential growth](../../img/sla_potential_growth.svg)
"""
function potential_growth(; SLA, nspecies, biomass, PAR)
    lais = calculate_LAI(; SLA, biomass)
    LAItot = sum(lais)

    if LAItot == 0
        return fill(0.0, nspecies), LAItot
    end

    RUE_max = 3 // 1000 # Maximum radiation use efficiency 3 g DM MJ-1
    α = 0.6   # Extinction coefficient, unitless

    total_growth = PAR * RUE_max * (1 - exp(-α * LAItot))
    species_growth = total_growth .* lais ./ LAItot

    return species_growth, LAItot
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
    return SLA .* biomass .* LAM
end

end # of module
