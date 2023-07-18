############################ Helper functions
function species_biomass(biomass)
    biomass_unit = unit(biomass[1])
    _, npatches, nspecies = size(biomass)
    mean_biomass = fill(0.0, nspecies) .* biomass_unit

    for pa in 1:npatches
        mean_patch_biomass = vec(
            mean(ustrip.(biomass[:, pa, :]), dims=1)
        ) .* biomass_unit

        mean_biomass += mean_patch_biomass
    end

    return mean_biomass
end

function total_t_biomass(biomass)
    biomass_unit = unit(biomass[1])
    return vec(sum(ustrip.(biomass[:, :, :]), dims=2:3)) .* biomass_unit
end

function patch_t_biomass(biomass)
    timesteps, npatches, _ = size(biomass)
    biomass_unit = unit(biomass[1])

    biomass_sum = sum(ustrip.(biomass[:, :, :]), dims=3) .* biomass_unit
    return reshape(biomass_sum, timesteps, npatches)
end

function test_date(x)
    return isnothing(
        tryparse(Dates.Date, x, Dates.dateformat"mm-dd")
    ) ? false : true
end
