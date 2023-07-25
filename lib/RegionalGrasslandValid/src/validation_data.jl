function get_validation_data(; plotID)
    soilmoisture_sub = @subset data.valid.soilmoisture :plotID.==plotID
    evaporation_sub = @subset data.valid.evaporation :plotID.==plotID
    satbiomass_sub = @subset data.valid.satbiomass :plotID.==plotID
    measuredveg_sub = @subset data.valid.measuredveg :plotID.==plotID

    return (;
        soilmoisture = soilmoisture_sub.soilmoisture,
        soilmoisture_t = soilmoisture_sub.sol_t,
        evaporation = evaporation_sub.evaporation,
        evaporation_t = evaporation_sub.sol_t,
        biomass = satbiomass_sub.biomass_kg_ha,
        biomass_t = satbiomass_sub.sol_t,
        measured_biomass = measuredveg_sub.biomass_kg_ha,
        measured_biomass1 = measuredveg_sub.biomass1_kg_ha,
        measured_biomass_t = measuredveg_sub.sol_t)
end
