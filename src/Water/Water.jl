module Water

using Unitful

"""
    change_water_reserve(;
        WR, precipitation, LAItot,
        PET, WHC, PWP)

TBW
"""
function change_water_reserve(;
    WR, precipitation, LAItot,
    PET, WHC, PWP)

    # Maximal measured value of PET
    # PETₘₐₓ= 8u"mm / d"

    # -------- Evapotranspiration
    AEv = evaporation(; WR, WHC, PET, LAItot)
    ATr = transpiration(; WR, PWP, WHC, PET, LAItot)
    AET = min(WR / u"d", ATr + AEv)

    # -------- Drainage
    drain = max(WR / u"d" + precipitation - WHC / u"d" - AET, 0u"mm / d")

    # -------- Total change in the water reserve
    water_change = precipitation - drain - AET

    return water_change
end


"""
    transpiration(;
    WR,
    PWP, WHC,
    PET,
    LAItot)

TBW
"""
function transpiration(;
    WR,
    PWP, WHC,
    PET,
    LAItot)

    # β₁ = 6.467 # unitless
    # β₂ = 7.623e-8 # unitless
    # exp_fun = -( β₂*PET/PETₘₐₓ + (1-PET/PETₘₐₓ)*β₁ )
    # @info W, (1 - exp(exp_fun*W) ) / (1-exp(exp_fun))
    # transpi = W * PET * min(1, LAI_tot/3) * (1 - exp(exp_fun*W) ) / (1-exp(exp_fun))

    # plant available water:
    W = max(0.0, (WR - PWP) / (WHC - PWP))

    transpi = W * PET * LAItot / 3

    return transpi
end

"""
    evaporation(; WR, WHC, PET, LAItot)

Aev(t) = WR(t) / WHC * PET(t)*[1 - min(1; LAI_tot(t)/3) ]
"""
function evaporation(; WR, WHC, PET, LAItot)
    W = WR / WHC
    return W * PET * (1 - min(1, LAItot / 3))
end

end # of module
