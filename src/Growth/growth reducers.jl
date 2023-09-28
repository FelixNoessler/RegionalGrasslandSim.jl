@doc raw"""
    radiation_reduction(; PAR, radiation_red)

Reduction of radiation use efficiency at light intensities higher than 5 ``MJ\cdot m^{-2}\cdot d^{-1}``

```math
\text{Rred} = \text{min}(1, 1-\gamma_1(\text{PAR}(t) - \gamma_2))
```
**comment to the equation/figure:** PAR values are usually between 0 and 15 ``MJ\cdot m^{-2}\cdot d^{-1}`` and therefore negative values of Rred are very unlikely
![Image of the radiation reducer function](../../img/radiation_reducer.svg)
"""
function radiation_reduction(; PAR, radiation_red)
    if !radiation_red
        @info "No radiation reduction!" maxlog=1
        return 1.0
    end

    γ1 = 0.0445u"m^2 * d / MJ" # Empirical parameter for a decrease
    # in RUE for high PAR values, m2 d MJ-1
    γ2 = 5.0u"MJ / (m^2 * d)" # Threshold value of PAR from which starts
    # a linear decrease in RUE, MJ m2 d-1

    return min(1.0, 1.0 − γ1 * (PAR − γ2))
end

"""
    temperature_reduction(; T, temperature_red)

TBW

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
        Waterred,
        sla_water,
        srsa_water,
        fun_response,
        WR,
        water_red,
        PET,
        PWP,
        WHC)

See for details: [Water stress](@ref water_stress)
"""
function water_reduction!(;
    Waterred,
    sla_water,
    srsa_water,
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

    ### ------------ species specific functional response
    sla_water_reduction!(; sla_water, fun_response, x)
    srsa_water_reduction!(; srsa_water, fun_response, x)

    @. Waterred = sla_water * srsa_water

    return nothing
end

"""
    sla_water_reduction!(;
        sla_water,
        fun_response,
        x)

TBW
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

TBW
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
        Nutred,
        amc_nut,
        srsa_nut,
        fun_response,
        nutrient_red,
        nutrients)

See for details: [Nutrient stress](@ref nut_stress)
"""
function nutrient_reduction!(;
    Nutred,
    amc_nut,
    srsa_nut,
    fun_response,
    nutrient_red,
    nutrients)
    if !nutrient_red
        @info "No nutrient reduction!" maxlog=1
        return 1.0
    end

    ### ------------ species specific functional response
    amc_nut_reduction!(; amc_nut, fun_response, x = nutrients)
    srsa_nut_reduction!(; srsa_nut, fun_response, x = nutrients)

    @. Nutred = max(amc_nut, srsa_nut)

    return nothing
end

"""
    amc_nut_reduction!(; amc_nut, fun_response, x)

TBW
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

TBW
"""
function srsa_nut_reduction!(; srsa_nut, fun_response, x)
    k_SRSA = 7

    @. srsa_nut = fun_response.srsa_nut_lower +
                  (fun_response.srsa_nut_upper - fun_response.srsa_nut_lower) /
                  (1 + exp(-k_SRSA * (x - fun_response.srsa_midpoint)))
    return nothing
end

"""
    seasonal_reduction(; ST, season_red)

TBW

![Image of the seasonal effect function](../../img/seasonal_reducer.svg)
"""
function seasonal_reduction(; ST, season_red)
    if !season_red
        @info "No seasonal reduction!" maxlog=1
        return 1.0
    end

    SEAₘᵢₙ = 0.67 # unitless
    SEAₘₐₓ = 1.33 # unitless

    ST₁ = 625  # u"°C d"
    ST₂ = 1300 # u"°C d"

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
