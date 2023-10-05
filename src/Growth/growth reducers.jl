@doc raw"""
    radiation_reduction(; PAR, radiation_red)

Reduction of radiation use efficiency at light intensities higher
than 5 ``MJ\cdot m^{-2}\cdot d^{-1}``

```math
\text{Rred} = \text{min}(1, 1-\gamma_1(\text{PAR}(t) - \gamma_2))
```

The equations and the parameter values are taken from [Schapendonk1998](@cite).

- `γ₁` is the empirical parameter for a decrease in RUE for high PAR values,
  here set to 0.0445 [m² d MJ⁻¹]
- `γ₂` is the threshold value of PAR from which starts a linear decrease in RUE,
  here set to 5 [MJ m⁻² d⁻¹]

comment to the equation/figure: PAR values are usually between 0 and
15 ``MJ\cdot m^{-2}\cdot d^{-1}`` and therefore negative values of
Rred are very unlikely
![Image of the radiation reducer function](../../img/radiation_reducer.svg)
"""
function radiation_reduction(; PAR, radiation_red)
    if !radiation_red
        @info "No radiation reduction!" maxlog=1
        return 1.0
    end

    γ1 = 0.0445u"m^2 * d / MJ"
    γ2 = 5.0u"MJ / (m^2 * d)"

    return min(1.0, 1.0 − γ1 * (PAR − γ2))
end

@doc raw"""
    temperature_reduction(; T, temperature_red)

Reduction of the potential growth if the temperature is low or too high
with a step function.

```math
\text{temperature_reduction}(T) =
    \begin{cases}
    0 & \text{if } T < T_0 \\
    \frac{T - T_0}{T_1 - T_0} & \text{if } T_0 < T < T_1 \\
    1 & \text{if } T_1 < T < T_2 \\
    \frac{T_3 - T}{T_3 - T_2} & \text{if } T_2 < T < T_3 \\
    0 & \text{if } T > T_3 \\
    \end{cases}
```

Equations are taken from [Moulin2021](@cite) and theses are based on
[Schapendonk1998](@cite). `T₁` is in [Moulin2021](@cite) a
species specific parameter, but here it is set to 12°C for all species.

- `T₀` is the lower temperature threshold for growth, here set to 3°C
- `T₁` is the lower bound for the optimal temperature for growth, here set to 12°C
- `T₂` is the upper bound for the optiomal temperature for growth, here set to 20°C
- `T₃` is the maximum temperature for growth, here set to 35°C

![Image of the temperature reducer function](../../img/temperature_reducer.svg)
"""
function temperature_reduction(; T, temperature_red)
    if !temperature_red
        @info "No temperature reduction!" maxlog=1
        return 1.0
    end

    T = ustrip(T)

    T₀ = 3  #u"°C"
    T₁ = 12 #u"°C"
    T₂ = 20 #u"°C"
    T₃ = 35 #u"°C"

    if T < T₀
        return 0
    elseif T < T₁
        return (T - T₀) / (T₁ - T₀)
    elseif T < T₂
        return 1
    elseif T < T₃
        return (T₃ - T) / (T₃ - T₂)
    else
        return 0
    end
end

"""
    water_reduction!(;
        calc,
        fun_response,
        WR,
        water_red,
        PET,
        PWP,
        WHC)

See for details: [Water stress](@ref water_stress)
"""
function water_reduction!(;
    calc,
    fun_response,
    WR,
    water_red,
    PET,
    PWP,
    WHC)
    if !water_red
        @info "No water reduction!" maxlog=1
        return 1.0
    end

    ## option 1: water reduction purely by water availability
    # x = WR > WHC ? 1.0 : WR > PWP ? (WR - PWP) / (WHC - PWP) : 0.0

    ## option 2: water reduction by water availability and
    ##           potential evapotranspiration
    W = WR > WHC ? 1.0 : WR > PWP ? (WR - PWP) / (WHC - PWP) : 0.0
    PETₘₐₓ = 8u"mm / d"
    β₁ = 6.467
    β₂ = 7.623e-8
    exp_fun = -(β₂ * PET / PETₘₐₓ + (1 - PET / PETₘₐₓ) * β₁)
    x = (1 - exp(exp_fun * W)) / (1 - exp(exp_fun))
    @. calc.water_splitted = x * calc.below_split
    ### ------------ species specific functional response
    sla_water_reduction!(;
        sla_water = calc.sla_water, fun_response, x = calc.water_splitted)
    srsa_water_reduction!(;
        srsa_water = calc.srsa_water, fun_response, x = calc.water_splitted)

    @. calc.Waterred = calc.sla_water * calc.srsa_water

    return nothing
end

"""
    sla_water_reduction!(;
        sla_water,
        fun_response,
        x)

Reduction of growth due to stronger water stress for higher specific leaf area (SLA).
"""
function sla_water_reduction!(;
    sla_water,
    fun_response,
    x)
    k_SLA = 5

    @. sla_water = fun_response.sla_water_lower +
                   (1 - fun_response.sla_water_lower) /
                   (1 + exp(-k_SLA * (x - fun_response.sla_water_midpoint)))

    return nothing
end

"""
    srsa_water_reduction!(; srsa_water, fun_response, x)

Reduction of growth due to stronger water stress for lower specific
root surface area per above ground biomass (`SRSA_above`).
"""
function srsa_water_reduction!(; srsa_water, fun_response, x)
    k_SRSA = 7

    @. srsa_water = fun_response.srsa_water_lower +
                    (fun_response.srsa_water_upper - fun_response.srsa_water_lower) /
                    (1 + exp(-k_SRSA * (x - fun_response.srsa_midpoint)))
    return nothing
end

"""
    nutrient_reduction!(;
        calc,
        fun_response,
        nutrient_red,
        nutrients)

See for details: [Nutrient stress](@ref nut_stress)
"""
function nutrient_reduction!(;
    calc,
    fun_response,
    nutrient_red,
    nutrients)
    if !nutrient_red
        @info "No nutrient reduction!" maxlog=1
        return 1.0
    end

    @. calc.nutrients_splitted = calc.below_split * nutrients

    ### ------------ species specific functional response
    amc_nut_reduction!(;
        amc_nut = calc.amc_nut, fun_response, x = calc.nutrients_splitted)
    srsa_nut_reduction!(;
        srsa_nut = calc.srsa_nut, fun_response, x = calc.nutrients_splitted)

    @. calc.Nutred = max(calc.amc_nut, calc.srsa_nut)

    return nothing
end

"""
    amc_nut_reduction!(; amc_nut, fun_response, x)

Reduction of growth due to stronger nutrient stress for lower
arbuscular mycorrhizal colonization (`AMC`).
"""
function amc_nut_reduction!(; amc_nut, fun_response, x)
    k_AMC = 7

    @. amc_nut = fun_response.myco_nut_lower +
                 (fun_response.myco_nut_upper - fun_response.myco_nut_lower) /
                 (1 + exp(-k_AMC * (x - fun_response.myco_nut_midpoint)))
    return nothing
end

"""
    srsa_nut_reduction!(; srsa_nut, fun_response, x)

Reduction of growth due to stronger nutrient stress for lower specific
root surface area per above ground biomass (`SRSA_above`).
"""
function srsa_nut_reduction!(; srsa_nut, fun_response, x)
    k_SRSA = 7

    @. srsa_nut = fun_response.srsa_nut_lower +
                  (fun_response.srsa_nut_upper - fun_response.srsa_nut_lower) /
                  (1 + exp(-k_SRSA * (x - fun_response.srsa_midpoint)))
    return nothing
end

@doc raw"""
    seasonal_reduction(; ST, season_red)

Reduction of growth due to seasonal effects. The function is based on
the yearly cumulative sum of the daily mean temperatures (`ST`).

```math
\text{seasonal}(ST) =
    \begin{cases}
    SEA_{min} & \text{if } ST < 200 \\
    SEAₘᵢₙ + (SEAₘₐₓ - SEAₘᵢₙ) * \frac{ST - 200}{ST₁ - 400} &
        \text{if } 200 < ST < ST₁ - 200 \\
    SEA_{max} & \text{if } ST₁ - 200 < ST < ST₁ - 100 \\
    SEAₘᵢₙ + (SEAₘᵢₙ - SEAₘₐₓ) * \frac{ST - ST₂}{ST₂ - ST₁ - 100} &
        \text{if } ST₁ - 100 < ST < ST₂ \\
    SEA_{min} & \text{if } ST > ST₂ \\
    \end{cases}
```

This empirical function was developed by [Jouven2006](@cite). In contrast to
[Jouven2006](@cite) `SEAₘᵢₙ`, `SEAₘₐₓ`, `ST₁` and `ST₂` are not species specific
parameters, but are fixed for all species. The values of the parameters are based on
[Jouven2006](@cite) and were chosen to resemble the mean of all functional
groups that were described there.

A seasonal factor greater than one means that growth is increased by the
use of already stored resources. A seasonal factor below one means that
growth is reduced as the plant stores resources [Jouven2006](@cite).

- `ST` is the yearly cumulative sum of the daily mean temperatures
- `SEAₘᵢₙ` is the minimum value of the seasonal effect, here set to 0.67 [-]
- `SEAₘₐₓ` is the maximum value of the seasonal effect, here set to 1.33 [-]
-  `ST₁` and `ST₂` are parameters that describe the thresholds of the step function,
   here set to 625 and 1300 [°C d]

![Image of the seasonal effect function](../../img/seasonal_reducer.svg)
"""
function seasonal_reduction(; ST, season_red)
    if !season_red
        @info "No seasonal reduction!" maxlog=1
        return 1.0
    end

    SEAₘᵢₙ = 0.7
    SEAₘₐₓ = 1.3
    ST₁ = 625
    ST₂ = 1300

    ST = ustrip(ST)

    if ST < 200
        return SEAₘᵢₙ
    elseif ST < ST₁ - 200
        return SEAₘᵢₙ + (SEAₘₐₓ - SEAₘᵢₙ) * (ST - 200) / (ST₁ - 400)
    elseif ST < ST₁ - 100
        return SEAₘₐₓ
    elseif ST < ST₂
        return SEAₘᵢₙ + (SEAₘᵢₙ - SEAₘₐₓ) * (ST - ST₂) / (ST₂ - (ST₁ - 100))
    else
        return SEAₘᵢₙ
    end
end
