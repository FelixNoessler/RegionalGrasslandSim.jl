function validation_input(;
    plotID,
    inf_p,
    nyears)
    nspecies = 25
    npatches = 1
    water_reduction = true
    nutrient_reduction = true

    initbiomass_sub = @subset data.input.initbiomass :plotID.==plotID
    clim_sub = @subset data.input.clim :plotID.==plotID
    pet_sub = @subset data.input.pet first.(:explo).==first(plotID)
    par_sub = @subset data.input.par first.(:explo).==first(plotID)
    nut_sub = @subset data.input.nut :plotID.==plotID
    soil_sub = @subset data.input.soil :plotID.==plotID
    mow_sub = @subset data.input.mow :plotID.==plotID
    graz_sub = @subset data.input.graz :plotID.==plotID

    ### ----------------- initial biomass
    initbiomass = initbiomass_sub.biomass_init[1] * u"kg / ha"

    ### ----------------- abiotic
    temperature = clim_sub.Ta_200 .* u"°C"
    temperature_sum = yearly_temp_cumsum(temperature)
    precipitation = clim_sub.precipitation .* u"mm / d"
    PAR = par_sub.PAR .* 10000 .* u"MJ / (d * ha)"
    PET = pet_sub.PET .* u"mm / d"
    nutrient_index = nut_sub.n_index[1]
    WHC = soil_sub.WHC[1] * u"mm"
    PWP = soil_sub.PWP[1] * u"mm"
    Clay = soil_sub.Clay[1]
    Silt = soil_sub.Silt[1]
    Sand = soil_sub.Sand[1]
    organic = soil_sub.organic[1]
    bulk = soil_sub.bulk[1]
    root_depth = soil_sub.root_depth[1]

    ### ----------------- traits
    traits = Traits.random_traits(nspecies;)
    sort!(traits, :SLA)
    relative_traits = Traits.relative_traits(; trait_data = traits)

    ### ----------------- mowing
    mowing_day = [mow_sub[year, "MowingDay$i"] for year in 1:nyears, i in 1:5]

    mowing_height = [mow_sub[year, "CutHeight_cm$i"] for year in 1:nyears, i in 1:5]
    height_missing = ismissing.(mowing_height) .&& .!ismissing.(mowing_day)

    mowing_height = convert(Matrix{Union{Missing, Int64}}, mowing_height)
    mowing_height[height_missing] .= 7
    height_should_be_missing = .!ismissing.(mowing_height) .&& ismissing.(mowing_day)
    mowing_height[height_should_be_missing] .= missing

    mowing_heights = [mowing_height[y, :][.!ismissing.(mowing_height[y, :])]
                      for y in 1:nyears]
    mowing_heights = convert(Vector{Vector{Int64}}, mowing_heights)

    mowing_days = [mowing_day[y, :][.!ismissing.(mowing_day[y, :])] for y in 1:nyears]
    mowing_days = convert(Vector{Vector{Int64}}, mowing_days)

    ### ----------------- grazing
    grazing_start = [graz_sub[year, "StartGrazingPeriod$i"] for year in 1:nyears, i in 1:4]
    days_grazing = [graz_sub[year, "DayGrazing$i"] for year in 1:nyears, i in 1:4]

    ### set the start day to missing, if grazing period == 0
    grazing_start = convert(Matrix{Union{Missing, Int64}}, grazing_start)
    grazing_start[iszero.(days_grazing)] .= missing

    ### calculate the end day of the grazing period
    grazing_end = grazing_start .+ days_grazing

    ### get the grazing intensity and transform 0 and NaNs to missing
    grazing_intensity = [graz_sub[year, "GrazingIntensity$i"] for year in 1:nyears, i in 1:4]
    grazing_intensity = convert(Matrix{Union{Missing, Float64}}, grazing_intensity)
    intensity_filter = isnan.(grazing_intensity) .|| iszero.(grazing_intensity)
    grazing_intensity[intensity_filter] .= missing

    ### convert matrices with missing to vectors of vectors
    grazing_start = [grazing_start[y, :][.!ismissing.(grazing_start[y, :])]
                     for y in 1:nyears]
    grazing_end = [grazing_end[y, :][.!ismissing.(grazing_end[y, :])] for y in 1:nyears]
    grazing_intensity = [grazing_intensity[y, :][.!ismissing.(grazing_intensity[y, :])]
                         for y in 1:nyears]

    ### derive the final grazing density vector
    grazing_density = fill(0.0u"ha ^ -1", 365 * nyears)
    for y in 1:nyears
        for i in eachindex(grazing_start[y])
            start_day = 365 * (y - 1) + grazing_start[y][i]
            end_day = 365 * (y - 1) + grazing_end[y][i]
            livestock_density = grazing_intensity[y][i]

            if isa(start_day, Number) && isa(end_day, Number) &&
               isa(livestock_density, Number)
                end_day = min(365 * nyears, end_day)

                grazing_density[start_day:end_day] .= livestock_density .* u"ha ^ -1"
            end
        end
    end

    return (inf_p,
        env_data = (;
            PAR,
            precipitation,
            temperature,
            temperature_sum,
            PET),
        site = (;
            initbiomass,
            nutrient_index,
            WHC,
            PWP,
            Clay,
            Silt,
            Sand,
            organic,
            bulk,
            root_depth),
        traits,
        relative_traits,
        mowing_days,
        mowing_heights,
        grazing = grazing_density,
        nspecies,
        npatches,
        water_reduction,
        nutrient_reduction,
        nyears)
end

function yearly_temp_cumsum(d)
    adj = 365
    nyears = length(d) ÷ adj + 1

    final_cumsum = Array{Float64}(undef, length(d))

    d = ustrip.(d)
    d[d .< 0.0] .= 0.0

    for y in 1:nyears
        sliced_d = d[(1 + adj * (y - 1)):min(adj * y, length(d))]
        final_cumsum[(1 + adj * (y - 1)):min(adj * y, length(d))] .= cumsum(sliced_d)
    end

    return final_cumsum
end
