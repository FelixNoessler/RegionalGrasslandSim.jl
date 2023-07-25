module Water

using Unitful

"""
    change_water_reserve()

TBW
"""
function change_water_reserve(;
    WR, precipitation, LAI_tot,
    PET, WHC, PWP)

    # Maximal measured value of PET
    # PETₘₐₓ= 8u"mm / d"

    # -------- Evapotranspiration
    AEv = evaporation(; WR, WHC, PET, LAI_tot)
    ATr = transpiration(; WR, PWP, WHC, PET, LAI_tot)
    AET = actual_evapotranspiration(; WR, ATr, AEv)

    # -------- Drainage
    drain = water_drainage(; WR, precipitation, WHC, AET)

    # -------- Total change in the water reserve
    water_change = precipitation - drain - AET

    return water_change, AEv
end

function actual_evapotranspiration(; WR, ATr, AEv)
    return min(WR * u"d^-1", ATr + AEv)
end

"""
    transpiration(;
        WR,
        PWP, WHC,
        PET,
        LAI_tot)

TBW
"""
function transpiration(;
    WR,
    PWP, WHC,
    PET,
    LAI_tot)

    # β₁ = 6.467 # unitless
    # β₂ = 7.623e-8 # unitless
    # exp_fun = -( β₂*PET/PETₘₐₓ + (1-PET/PETₘₐₓ)*β₁ )
    # @info W, (1 - exp(exp_fun*W) ) / (1-exp(exp_fun))
    # transpi = W * PET * min(1, LAI_tot/3) * (1 - exp(exp_fun*W) ) / (1-exp(exp_fun))

    # plant available water:
    W = max(0.0, (WR - PWP) / (WHC - PWP))

    transpi = W * PET * LAI_tot / 3

    return transpi
end

"""
    evaporation(; WR, WHC, PET, LAI_tot)

Aev(t) = WR(t) / WHC * PET(t)*[1 - min(1; LAI_tot(t)/3) ]
"""
function evaporation(; WR, WHC, PET, LAI_tot)
    W = WR / WHC
    return W * PET * (1 - min(1, LAI_tot / 3))
end

"""
    water_drainage(; WR, P, WHC, AET)

Δ(𝑡) = max(𝑊𝑅(𝑡) + 𝑃 (𝑡) − 𝐴𝐸𝑇 (𝑡) − 𝑊𝐻𝐶 ; 0)
"""
function water_drainage(; WR, precipitation, WHC, AET)
    return max(WR * u"d^-1" + precipitation - WHC * u"d^-1" - AET,
        0.0u"mm / d")
end

end # of module
