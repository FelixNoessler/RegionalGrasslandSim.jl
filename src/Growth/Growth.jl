module Growth

using Unitful

include("growth reducers.jl")
include("defoliation.jl")
include("senescence.jl")

"""
    growth(; t, p, calc, biomass, WR)

Calculates the actual growth of the plant species.
"""
function growth(; t, p, calc, biomass, WR)
    #### potential growth
    calc.LAItot = potential_growth(;
        pot_growth = calc.pot_growth,
        LAIs = calc.LAIs,
        SLA = p.species.SLA,
        nspecies = p.nspecies,
        biomass,
        PAR = p.env_data.PAR[t])

    ### influence of the height of plants
    calc.CHinfluence .= height_influence(; biomass, CH = p.species.CH)

    #### below ground competition --> trait similarity and abundance
    calc.below = below_ground_competition(;
        biomass,
        trait_similarity = p.trait_similarity,
        nspecies = p.nspecies,
        below_competition_strength = p.inf_p.below_competition_strength)

    #### growth reducer
    calc.Waterred .= water_reduction(;
        WR,
        PWP = p.site.PWP,
        WHC = p.site.WHC,
        fun_response = p.species.fun_response,
        water_red = p.water_reduction,
        PET = p.env_data.PET[t],
        max_SRSA_water_reduction = p.inf_p.max_SRSA_water_reduction,
        max_SLA_water_reduction = p.inf_p.max_SLA_water_reduction)
    calc.Nutred .= nutrient_reduction(;
        fun_response = p.species.fun_response,
        nutrient_red = p.nutrient_reduction,
        nutrients = p.site.nutrient_index,
        max_AMC_nut_reduction = p.inf_p.max_AMC_nut_reduction,
        max_SRSA_nut_reduction = p.inf_p.max_SRSA_nut_reduction)
    Rred = radiation_reduction(; PAR = p.env_data.PAR[t])
    Tred = temperature_reduction(; T = p.env_data.temperature[t])
    Seasonalred = seasonal_reduction(; ST = p.env_data.temperature_sum[t])


    @. calc.species_specific_red = calc.CHinfluence * calc.below * calc.Waterred * calc.Nutred
    reduction = Rred * Tred  * Seasonalred

    #### final growth
    calc.act_growth .= calc.pot_growth .* reduction .* calc.species_specific_red
    if any(calc.act_growth .< 0u"kg / (ha * d)")
        @error "act_growth below zero: $(calc.act_growth)" maxlog=10
    end

    return nothing
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
    biomass = ustrip.(biomass)

    for i in 1:nspecies
        trait_sim = @view trait_similarity[i, :]
        x = sum(trait_sim .* biomass)
        reduction_coefficient[i] = exp(-below_competition_strength / 1000 * x)
    end

    return reduction_coefficient
end

@doc raw"""
    potential_growth(; pot_growth, LAIs, SLA, nspecies, biomass, PAR)

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
function potential_growth(; pot_growth, LAIs, SLA, nspecies, biomass, PAR)
    LAItot = calculate_LAI(; LAIs, SLA, biomass)
    if LAItot < 0
        @error "LAItot below zero: $LAItot" maxlog=10
    end

    if LAItot == 0
        return fill(0.0, nspecies)u"kg / ha / d"
    end

    RUE_max = 3 // 1000 * u"kg / MJ" # Maximum radiation use efficiency 3 g DM MJ-1
    α = 0.6   # Extinction coefficient, unitless
    pot_growth_tot = PAR * RUE_max * (1 - exp(-α * LAItot))
    @. pot_growth = pot_growth_tot * LAIs / LAItot

    return LAItot
end

@doc raw"""
    community_weighted_mean_height(; biomass, CH)

```math
\text{CH}_{\text{cwm}} =
    \frac{\sum \text{biomass}
    \cdot \text{CH}}{\sum \text{biomass}}
```
"""
function community_weighted_mean_height(; biomass, CH)
    return sum(biomass .* CH) / sum(biomass)
end

@doc raw"""
    height_influence(; biomass, CH, CH_strength = 0.5)

```math
\text{CHinfluence} =
    1 +
    \frac{\text{CH}\cdot\text{CH}_{\text{strength}}}{\text{CH}_{\text{cwm}}}
    -\text{CH}_{\text{strength}}
```

- `CH_strength` lies between 0 (no influence) and 1 (strong influence of the plant height)
- the community weighted mean height `CH_cwm` is calculated by [`community_weighted_mean_height`](@ref)

In these plots all three plant species have an equal biomass:
![](../../img/ch_influence_05.svg)
![](../../img/ch_influence_08.svg)
"""
function height_influence(; biomass, CH, CH_strength = 0.5)
    CH_cwm = community_weighted_mean_height(; biomass, CH)
    return @. CH * CH_strength / CH_cwm - CH_strength + 1
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
    LAIs .= uconvert.(NoUnits, SLA .* biomass * LAM)
    return sum(LAIs)
end

end # of module
