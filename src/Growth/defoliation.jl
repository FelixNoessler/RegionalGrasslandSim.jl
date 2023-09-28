@doc raw"""
    mowing!(;
        calc,
        mowing_height,
        days_since_last_mowing,
        height,
        biomass,
        mowing_mid_days)

```math
\begin{align}
    \lambda &= \frac{\text{mown_height}}{height}\\
    \text{mow_factor} &= \frac{1}{1+exp(-0.1*(\text{days_since_last_mowing}-\text{mowing_mid_days})}\\
    \text{mow} &= \lambda \cdot \text{biomass}
\end{align}
```

The mow_factor has been included to account for the fact that less biomass is mown
when the last mowing event was not long ago.
Influence of mowing for plant species with different heights ($height$):
![Image of mowing effect](../../img/mowing.svg)

Visualisation of the `mow_factor`:
![](../../img/mow_factor.svg)
"""
function mowing!(;
    calc,
    mowing_height,
    days_since_last_mowing,
    height,
    biomass,
    mowing_mid_days)

    # --------- mowing parameter λ
    calc.mown_height .= height .- mowing_height
    calc.mown_height .= max.(calc.mown_height, 0.01u"m")
    calc.mowing_λ .= calc.mown_height ./ height

    # --------- if meadow is too often mown, less biomass is removed
    ## the 'mowing_mid_days' is the day where the plants are grown
    ## back to their normal size/2
    mow_factor = 1 / (1 + exp(-0.05 * (days_since_last_mowing - mowing_mid_days)))

    # --------- add the removed biomass to the defoliation vector
    calc.defoliation .+= mow_factor .* calc.mowing_λ .* biomass .* u"d^-1"

    return nothing
end

@doc raw"""
    grazing!(; calc, LD, biomass, ρ, grazing_half_factor)

```math
\begin{align}
μₘₐₓ &= κ \cdot \text{LD} \\
h &= \frac{1}{μₘₐₓ} \\
a &= \frac{1}{\text{grazing_half_factor}^2 \cdot h} \\
\text{totgraz} &= \frac{a \cdot (\sum \text{biomass})^2}{1 + a\cdot h\cdot (\sum \text{biomass})^2} \\
\text{graz} &= \text{totgraz} \cdot \frac{ρ \cdot \text{biomass}}{\sum ρ \cdot \text{biomass}}
\end{align}
```

- `LD` daily livestock density (livestock units ha⁻¹)
- `κ` daily consumption of one livestock unit, follows [Gillet2008](@cite)
- `ρ` appetence of the plant species for livestock, dependent on nitrogen per leaf mass (LNCM)
- `grazing_half_factor` is the half-saturation constant
- equation partly based on [Moulin2021](@cite)

Influence of grazing (livestock density = 2), all plant species have an equal amount of biomass (total biomass / 3):
![Image of grazing effect](../../img/grazing.svg)

Influence of `grazing_half_factor` (`LD` is set to 2):
![](../../img/grazing_half_factor.svg)
"""
function grazing!(; calc, LD, biomass, ρ, grazing_half_factor)
    κ = 22u"kg / d"
    k_exp = 2
    μₘₐₓ = κ * LD
    h = 1 / μₘₐₓ
    a = 1 / (grazing_half_factor^k_exp * h)

    ## Exponentiation of Quantity with a variable is type unstable
    ## therefore this is a workaround, k_exp = 2
    # https://painterqubits.github.io/Unitful.jl/stable/trouble/#Exponentiation
    biomass_exp = sum(biomass)^2

    total_grazed_biomass = a * biomass_exp / (1u"kg / ha"^k_exp + a * h * biomass_exp)

    @. calc.biomass_ρ = ρ * biomass
    calc.grazed_share .= calc.biomass_ρ ./ sum(calc.biomass_ρ)

    #### add grazed biomass to defoliation
    @. calc.defoliation += calc.grazed_share * total_grazed_biomass

    return nothing
end

@doc raw"""
    trampling!(; calc, LD, biomass, height, trampling_factor)

```math
\begin{align}
ω &= \frac{\text{trampling_factor}}{height^{0.25}} \\
\text{trampled_biomass} &=
    \begin{cases}
        \text{biomass},
            & \text{if LD > ω}\\
        \frac{\text{biomass}}{2}
        \cdot \left(1- cos\left(\frac{π\cdot\text{LD}}{ω}\right)\right),
            & \text{otherwise}
	\end{cases}
\end{align}
```

It is assumed that tall plants (trait: $height$) are stronger affected by trampling.
A cosine function is used to model the influence of trampling.

If the livestock density is higher than $ω$, all the biomass of that plant
species will be removed. This is unlikely to be the case.

- biomass [$\frac{kg}{ha}$]
- LD daily livestock density [$\frac{\text{livestock units}}{ha}$]
- trampling_factor [$ha$]
- height canopy height [$m$]

![Image of trampling effect](../../img/trampling.svg)
"""
function trampling!(; calc, LD, biomass, height, trampling_factor)
    ## higher values of the trampling factor with unit: "m^0.25 / ha"
    ## --> less loss due to trampling
    ## values shouldn't be lower than 50,
    ## otherwise (larger) plants are too heavily influenced
    @. calc.trampling_ω = trampling_factor / ustrip(height)^0.25 * u"1/ha"
    @. calc.trampled_biomass = biomass * 0.5 * (1 - cos(π * LD / calc.trampling_ω))

    @. calc.trampling_high_LD = LD > calc.trampling_ω
    if any(calc.trampling_high_LD)
        @warn """
              trampling removed all biomass of at least one plant species
              during one day:
              - lifestock density LD=$(round(ustrip(LD); digits=2))
              - trampling_factor=$(round(ustrip(trampling_factor); digits=2))
              """ maxlog=3
        calc.trampled_biomass[calc.trampling_high_LD] .= biomass[calc.trampling_high_LD]
    end
    calc.defoliation .+= calc.trampled_biomass ./ u"d"
    return nothing
end
