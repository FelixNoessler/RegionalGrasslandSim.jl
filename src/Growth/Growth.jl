module Growth

using Unitful
using Distributions
import LinearAlgebra

include("growth reducers.jl")
include("defoliation.jl")
include("senescence.jl")

"""
    growth!(; t, p, calc, biomass, WR)

Calculates the actual growth of the plant species.
"""
function growth!(; t, p, calc, biomass, WR, nutrients, WHC, PWP)
    #### potential growth
    LAItot = potential_growth!(;
        calc,
        potgrowth_included = p.included.potgrowth_included,
        SLA = p.species.SLA,
        biomass,
        PAR = p.daily_data.PAR[t])

    ### influence of the height of plants
    height_influence!(;
        calc,
        biomass,
        height_included = p.included.height_included,
        height = p.species.height)

    #### below ground competition --> trait similarity and abundance
    below_ground_competition!(;
        calc,
        biomass,
        below_included = p.included.below_included,
        trait_similarity = p.trait_similarity,
        belowground_density_effect = p.inf_p.belowground_density_effect)

    #### growth reducer
    water_reduction!(;
        calc,
        WR,
        PWP,
        WHC,
        fun_response = p.species.fun_response,
        water_red = p.included.water_red,
        PET = p.daily_data.PET[t])
    nutrient_reduction!(;
        calc,
        fun_response = p.species.fun_response,
        nutrient_red = p.included.nutrient_red,
        nutrients)
    Rred = radiation_reduction(;
        PAR = p.daily_data.PAR[t],
        radiation_red = p.included.radiation_red)
    Tred = temperature_reduction(;
        T = p.daily_data.temperature[t],
        temperature_red = p.included.temperature_red)
    Seasonalred = seasonal_reduction(;
        ST = p.daily_data.temperature_sum[t],
        season_red = p.included.season_red)

    calc.species_specific_red .= calc.heightinfluence .* calc.Waterred .*
                                 calc.Nutred #.* calc.below
    reduction = Rred * Tred * Seasonalred

    #### final growth
    @. calc.act_growth .= calc.pot_growth * reduction * calc.species_specific_red

    @. calc.neg_act_growth = calc.act_growth < 0u"kg / (ha * d)"
    if any(calc.neg_act_growth)
        @error "act_growth below zero: $(calc.act_growth)" maxlog=10
    end

    return LAItot
end

@doc raw"""
    below_ground_competition!(;
        below,
        traitsimilarity_biomass,
        biomass, below_included,
        trait_similarity,
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

The `trait_similarity` is computed before the start of the simulation
([calculation of traits similarity](@ref trait_similarity)).
and includes the traits arbuscular mycorrhizal colonisation rate (`AMC`)
and the root surface area devided by the above ground biomass (`SRSA_above`).
"""
function below_ground_competition!(;
    calc,
    biomass, below_included,
    trait_similarity,
    belowground_density_effect)
    if !below_included
        @info "No below ground competition!" maxlog=1
        @. calc.below_split = 1.0
        return nothing
    end

    LinearAlgebra.mul!(calc.traitsimilarity_biomass, trait_similarity, biomass)

    x = calc.traitsimilarity_biomass
    a = 7.5
    b = (log(10.0^belowground_density_effect) - a) / 1000
    @. calc.below_split = exp(a + b * x) / (1 + exp(a + b * x))

    return nothing
end

@doc raw"""
    potential_growth!(; calc, SLA, biomass, PAR, potgrowth_included)

Calculates the potential growth of all plant species
in a specific patch.

This function is called each time step (day) for each patch.
The `PAR` value is the photosynthetically
active radiation of the day.

First, the leaf area indices of all species are calculated
(see [`Growth.calculate_LAI`](@ref)). Then, the total leaf area is
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
function potential_growth!(; calc, SLA, biomass, PAR, potgrowth_included)
    LAItot = calculate_LAI(; LAIs = calc.LAIs, SLA, biomass)
    if LAItot < 0
        @error "LAItot below zero: $LAItot" maxlog=10
    end

    if LAItot == 0 || !potgrowth_included
        @info "Zero potential growth!" maxlog=1
        calc.pot_growth .= 0.0u"kg / (ha * d)"
        return LAItot
    end

    RUE_max = 3 // 1000 * u"kg / MJ" # Maximum radiation use efficiency 3 g DM MJ-1
    α = 0.6   # Extinction coefficient, unitless
    pot_growth_tot = PAR * RUE_max * (1 - exp(-α * LAItot))
    @. calc.pot_growth = pot_growth_tot * calc.LAIs / LAItot

    return LAItot
end

@doc raw"""
    community_weighted_mean_height(; calc, biomass, height)

```math
\text{height}_{\text{cwm}} =
    \frac{\sum \text{biomass}
    \cdot \text{height}}{\sum \text{biomass}}
```
"""
function community_weighted_mean_height(; calc, biomass, height)
    calc.biomass_height .= biomass .* height
    d = sum(calc.biomass_height) / sum(biomass)
    return d
end

@doc raw"""
    height_influence!(;
        calc, biomass, height, height_included, height_strength = 0.5)

```math
\text{heightinfluence} =
    1 +
    \frac{\text{height}\cdot\text{height}_{\text{strength}}}{\text{height}_{\text{cwm}}}
    -\text{height}_{\text{strength}}
```

- `height_strength` lies between 0 (no influence) and 1 (strong influence of the plant height)
- the community weighted mean height `height_cwm` is calculated
  by [`Growth.community_weighted_mean_height`](@ref)

In these plots all three plant species have an equal biomass:
![](../../img/height_influence_01.svg)
![](../../img/height_influence_05.svg)
"""
function height_influence!(;
    calc, biomass, height, height_included, height_strength)
    if !height_included
        @info "Height influence turned off!" maxlog=1
        return 1.0
    end

    height_cwm = community_weighted_mean_height(; calc, biomass, height)
    @. calc.heightinfluence = height * height_strength / height_cwm - height_strength + 1.0

    return nothing
end

@doc raw"""
    calculate_LAI(; SLA, biomass, LAIs)

Calculate the leaf area index of all species of one habitat patch.

```math
\begin{align}
\text{LAI} &= \text{SLA} \cdot \text{biomass} \cdot \text{LAM} \\
\text{LAI}_{\text{tot}} &= \sum \text{LAI}
\end{align}
```

- `SLA` specific leaf area [m² g⁻¹]
- `LAM` Proportion of laminae in green biomass [unitless], the value 0.62 is derived by [Jouven2006](@cite)
- `biomass` [kg ha⁻¹]

There is a unit conversion from the `SLA` and the `biomass`
to the unitless `LAI` involved.

The array `LAIs` is mutated inplace.
"""
function calculate_LAI(; SLA, biomass, LAIs)
    LAM = 0.62 # Proportion of laminae in green biomass
    LAIs .= uconvert.(NoUnits, SLA .* biomass .* LAM)
    return sum(LAIs)
end

end # of module
