function calc_CWM(; biomass, trait_data)
    tend, npatches, nspecies  = size(biomass)
    trait_data = ustrip.(trait_data)
    biomass = ustrip.(biomass)
    cwm = Array{Float64}(undef, tend)

    for t in 1:tend
        for p in 1:npatches
            biomass_t_p = @view biomass[t, p, :]
            total_biomass = sum(biomass_t_p)
            weights = biomass_t_p ./ total_biomass
            cwm[t] = sum(weights .* trait_data)
        end
    end

    return cwm
end


function calc_CWV(; biomass, trait_data)
    tend, npatches, nspecies  = size(biomass)
    trait_data = ustrip.(trait_data)
    biomass = ustrip.(biomass)

    mean_trait = mean(trait_data)
    cwv = Array{Float64}(undef, tend)

    for t in 1:tend
        for p in 1:npatches
            biomass_t_p = @view biomass[t, p, :]
            total_biomass = sum(biomass_t_p)
            weights = biomass_t_p ./ total_biomass
            cwv[t] = sum(weights .* (trait_data .- mean_trait) .^ 2)
        end
    end

    return cwv
end
