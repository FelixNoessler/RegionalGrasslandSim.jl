@doc raw"""
    radiation_reduction(; PAR)

Reduction of radiation use efficiency at light intensities higher than 5 ``MJ\cdot m^{-2}\cdot d^{-1}``

```math
\text{Rred} = \text{min}(1, 1-\gamma_1(\text{PAR}(t) - \gamma_2))
```
**comment to the equation/figure:** PAR values are usually between 0 and 15 ``MJ\cdot m^{-2}\cdot d^{-1}`` and therefore negative values of Rred are very unlikely
![Image of the radiation reducer function](../../img/radiation_reducer.svg)
"""
function radiation_reduction(; PAR)
    γ1 = 0.0445u"m^2 * d / MJ" # Empirical parameter for a decrease
    # in RUE for high PAR values, m2 d MJ-1
    γ2 = 5.0u"MJ / (m^2 * d)" # Threshold value of PAR from which starts
    # a linear decrease in RUE, MJ m2 d-1

    return min(1.0, 1.0 − γ1 * (PAR − γ2))
end

"""
    temperature_reduction(; T)

TBW

![Image of the temperature reducer function](../../img/temperature_reducer.svg)
"""
function temperature_reduction(; T)
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
    water_reduction(;
        fun_response,
        WR,
        water_red,
        max_SLA_water_reduction,
        max_SRSA_water_reduction,
        PET,
        PWP,
        WHC)

See for details: [Water stress](@ref water_stress)
"""
function water_reduction(;
    fun_response,
    WR,
    water_red,
    max_SLA_water_reduction,
    max_SRSA_water_reduction,
    PET,
    PWP,
    WHC)

    if !water_red
        return 1.0
    end

    x = WR > WHC ? 1.0 : WR > PWP ? (WR - PWP) / (WHC - PWP) : 0.0

    ## could be added: on days with high PET, more water is needed
    # PETₘₐₓ= 8u"mm / d",
    # β₁=6.467,
    # β₂=7.623e-8
    # exp_fun = -(β₂*PET/PETₘₐₓ + (1-PET/PETₘₐₓ)*β₁)
    # x = (1 - exp(exp_fun*W) ) / (1-exp(exp_fun))

    ### ------------ species specific functional response
    sla_y = sla_water_reduction(; fun_response, x, max_SLA_water_reduction)
    srsa_y = srsa_water_reduction(; fun_response, x, max_SRSA_water_reduction)

    return sla_y .* srsa_y
end

"""
    sla_water_reduction(; fun_response, x, max_SLA_water_reduction)

TBW
"""
function sla_water_reduction(;
    fun_response,
    x,
    max_SLA_water_reduction)

    k_SLA = 5
    A_SLA = 1-max_SLA_water_reduction
    return @. A_SLA + (1 - A_SLA) /
            (1 + exp(-k_SLA * (x - fun_response.sla_water_midpoint)))
end

"""
    srsa_water_reduction(; fun_response, x, max_SRSA_water_reduction)

TBW
"""
function srsa_water_reduction(; fun_response, x, max_SRSA_water_reduction)
    k_SRSA = 7
    A_SRSA = 1-max_SRSA_water_reduction
    K_SRSA_prep = fun_response.srsa_right_bound
    K_SRSA = @. K_SRSA_prep + (1 - K_SRSA_prep) * A_SRSA

    return @. A_SRSA +  (K_SRSA  - A_SRSA) /
            (1 + exp(-k_SRSA * (x - fun_response.srsa_midpoint)))
end


"""
    nutrient_reduction(;
        fun_response,
        nutrient_red,
        nutrients,
        max_AMC_nut_reduction,
        max_SRSA_nut_reduction)

See for details: [Nutrient stress](@ref nut_stress)
"""
function nutrient_reduction(;
    fun_response,
    nutrient_red,
    nutrients,
    max_AMC_nut_reduction,
    max_SRSA_nut_reduction)

    if !nutrient_red
        return 1.0
    end

    ### ------------ species specific functional response
    amc_red = amc_nut_reduction(; fun_response, x=nutrients, max_AMC_nut_reduction)
    srsa_red = srsa_nut_reduction(; fun_response, x=nutrients, max_SRSA_nut_reduction)

    return max.(amc_red, srsa_red)
end


"""
    amc_nut_reduction(; fun_response, x, max_AMC_nut_reduction)

TBW
"""
function amc_nut_reduction(; fun_response, x, max_AMC_nut_reduction)
    k_AMC = 7
    A_AMC = 1-max_AMC_nut_reduction
    K_AMC_prep = fun_response.myco_nutr_right_bound
    K_AMC = @. K_AMC_prep + (1 - K_AMC_prep) * A_AMC

    return @. A_AMC + (K_AMC - A_AMC) /
            (1 + exp(-k_AMC * (x - fun_response.myco_nutr_midpoint)))
end

"""
    srsa_nut_reduction(; fun_response, x, max_SRSA_nut_reduction)

TBW
"""
function srsa_nut_reduction(; fun_response, x, max_SRSA_nut_reduction)
    k_SRSA = 7
    A_SRSA = 1-max_SRSA_nut_reduction
    K_SRSA_prep = fun_response.srsa_right_bound
    K_SRSA = @. K_SRSA_prep + (1 - K_SRSA_prep) * A_SRSA

    return @. A_SRSA + (K_SRSA - A_SRSA) /
            (1 + exp(-k_SRSA * (x - fun_response.srsa_midpoint)))
end



"""
    seasonal_reduction()

TBW

![Image of the seasonal effect function](../../img/seasonal_reducer.svg)
"""
function seasonal_reduction(;
        ST # accumulated degree days
    )
    SEAₘᵢₙ = 0.67 # unitless
    SEAₘₐₓ = 1.33 # unitless

    ST₁ = 625  # u"°C d"
    ST₂ = 1300 # u"°C d"

    if ST < 200
        return SEAₘᵢₙ
    elseif ST < ST₁-200
        return SEAₘᵢₙ + (SEAₘₐₓ-SEAₘᵢₙ)*(ST-200)/(ST₁-400)
    elseif ST < ST₁-100
        return SEAₘₐₓ
    elseif ST < ST₂
        return SEAₘᵢₙ + (SEAₘᵢₙ-SEAₘₐₓ)*(ST-ST₂)/(ST₂-(ST₁-100))
    else
        return SEAₘᵢₙ
    end
end
