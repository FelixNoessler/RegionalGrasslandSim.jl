@doc raw"""
    mowing(;
        mowing_height,
        days_since_last_mowing,
        CH,
        biomass,
        mowing_mid_days)

```math
\begin{align}
    \lambda &= \frac{\text{mown_height}}{CH}\\
    \text{mow_factor} &= \frac{1}{1+exp(-0.1*(\text{days_since_last_mowing}-\text{mowing_mid_days})}\\
    \text{mow} &= \lambda \cdot \text{biomass}
\end{align}
```

The mow_factor has been included to account for the fact that less biomass is mown
when the last mowing event was not long ago.
Influence of mowing for plant species with different heights ($CH$):
![Image of mowing effect](../../img/mowing.svg)

Visualisation of the `mow_factor`:
![](../../img/mow_factor.svg)
"""
function mowing(;
    mowing_height,
    days_since_last_mowing,
    CH,
    biomass,
    mowing_mid_days)

    # --------- mowing parameter λ
    mown_height = CH .- mowing_height / 100 * u"m"
    mown_height = max.(mown_height, 0.01u"m")
    λ = mown_height ./ CH

    # --------- if meadow is too often mown, less biomass is removed
    ## the 'mowing_mid_days' is the day where the plants are grown
    ## back to their normal size/2
    mow_factor = 1/(1+exp(-0.05*(days_since_last_mowing-mowing_mid_days)))

    # --------- biomass that is removed by mowing
    removed_biomass = mow_factor .* λ .* biomass .* u"d^-1"

    return removed_biomass
end


@doc raw"""
    grazing(; LD, biomass, ρ, nspecies, grazing_half_factor)

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
function grazing(; LD, biomass, ρ, nspecies, grazing_half_factor)
    if iszero(LD)
        return zeros(nspecies)u"kg / (ha * d)"
    end

    κ = 22u"kg / d"

    k_exp = 2
    μₘₐₓ = κ * LD
    h = 1 / μₘₐₓ
    a = 1 / (grazing_half_factor^k_exp * h)

    graz = a * sum(biomass)^k_exp  / (1u"kg / ha" ^ k_exp + a*h*sum(biomass)^k_exp)
    share = (ρ .* biomass) ./ sum(ρ .* biomass)

    return graz .* share
end


@doc raw"""
    trampling(; LD, biomass, LA, CH, trampling_factor)

```math
\begin{align}
ω &= \frac{\text{trampling_factor}}{CH^{0.25}} \\
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

It is assumed that tall plants (trait: $CH$) are stronger affected by trampling.
A cosine function is used to model the influence of trampling.

If the livestock density is higher than $ω$, all the biomass of that plant
species will be removed. This is unlikely to be the case.

- biomass [$\frac{kg}{ha}$]
- LD daily livestock density [$\frac{\text{livestock units}}{ha}$]
- trampling_factor [$ha$]
- CH canopy height [$m$]

![Image of trampling effect](../../img/trampling.svg)
"""
function trampling(; LD, biomass, CH, nspecies, trampling_factor)
    if iszero(LD)
        return zeros(nspecies)u"kg/(ha*d)"
    end

    ## higher values of the trampling factor
    ## --> less loss due to trampling
    ## values shouldn't be lower than 50,
    ## otherwise (larger) plants are too heavily influenced
    ω = @. trampling_factor / ustrip(CH)^0.25 * u"1/ha"
    trampled_biomass = @. biomass * 0.5 * (1 - cos(π*LD/ω))

    if any(LD .> ω)
        @warn "Very high lifestock densities (LD=$LD)!"
        high_LD = LD .> ω
        trampled_biomass[high_LD] .= biomass[high_LD]
    end

    return trampled_biomass .* u"d^-1"
end
