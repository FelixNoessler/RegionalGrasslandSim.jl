module Growth


"""
    greet_your_package_name()

Hello!
"""
function greet_your_package_name()
    return "Hello RegionalGrasslandSim!"
end

function potential_growth(p, biomass; PAR)
    lais = calculate_LAI(p, biomass)
    LAItot = sum(lais)
    
 
    RUE_max = 3
    α = 0.6

    total_growth = 10 * 0.0001 * PAR * RUE_max * (1 -  exp(-α * LAItot))

    species_growth = total_growth .* lais ./ LAItot

    return species_growth
end

function calculate_LAI(p, biomass)
    LAM = 0.62
    return p.species.sla .* biomass .* LAM .* 0.1
end

function senescence(; biomass, t, p)
    # include a seasonal effect
    # less senescence in spring, 
    # high senescens rate in autumn 
    # Ψ₁, Ψ₂ = 775, 3000
    # SEN_min, SEN_max = 1, 3

    return 0.01 .* p.species.μ .* biomass
end



end # of module