module Growth


"""
    greet_your_package_name()

Hello!
"""
function greet_your_package_name()
    return "Hello RegionalGrasslandSim!"
end

function potential_growth(p, biomass)
    lais = calculate_LAI(p, biomass)
    LAItot = calculate_totLAI(lais)
    
    PAR = 10
    RUE_max = 3
    α = 0.6

    total_growth = 10 * PAR * RUE_max * (1 -  exp(-α * LAItot))

    species_growth = total_growth .* lais ./ LAItot

    return species_growth
end

function calculate_LAI(p, biomass)
    LAM = 0.62
    return p.species.sla .* biomass  .* LAM .* 0.1
end

function calculate_totLAI(lai_values)
    return sum(lai_values)
end

end # of module