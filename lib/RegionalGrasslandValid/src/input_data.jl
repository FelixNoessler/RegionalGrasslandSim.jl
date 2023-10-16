function validation_input_plots(;
    plotIDs,
    nspecies,
    startyear, endyear,
    npatches = 1,
    nutheterog,
    senescence_included = true,
    potgrowth_included = true,
    mowing_included = true,
    grazing_included = true,
    below_included = true,
    height_included = true,
    water_red = true,
    nutrient_red = true,
    temperature_red = true,
    season_red = true,
    radiation_red = true,
    constant_init_biomass = true,
    constant_seed = false)
    static_inputs = Dict()

    for plotID in plotIDs
        static_inputs[plotID] = validation_input(;
            plotID,
            nspecies,
            startyear, endyear,
            npatches,
            nutheterog,
            senescence_included,
            potgrowth_included,
            mowing_included,
            grazing_included,
            below_included,
            height_included,
            water_red,
            nutrient_red,
            temperature_red,
            season_red,
            radiation_red,
            constant_init_biomass,
            constant_seed)
    end

    return static_inputs
end

function validation_input(;
    plotID,
    nspecies,
    startyear, endyear,
    nutheterog,
    npatches = 1,
    senescence_included = true,
    potgrowth_included = true,
    mowing_included = true,
    grazing_included = true,
    below_included = true,
    height_included = true,
    water_red = true,
    nutrient_red = true,
    temperature_red = true,
    season_red = true,
    radiation_red = true,
    constant_init_biomass = true,
    constant_seed = false)
    yearrange = Ref(startyear:endyear)

    ###### daily data
    clim_sub = @subset data.input.clim :plotID .== plotID.&&:year .∈ yearrange
    pet_sub = @subset data.input.pet first.(:explo) .== first(plotID).&&:year .∈ yearrange
    par_sub = @subset data.input.par first.(:explo) .== first(plotID).&&:year .∈ yearrange
    mow_sub = @subset data.input.mow :plotID .== plotID.&&:year .∈ yearrange
    graz_sub = @subset data.input.graz :plotID .== plotID.&&:year .∈ yearrange

    daily_data_prep = @chain clim_sub begin
        innerjoin(_, pet_sub, on = :date, makeunique = true)
        innerjoin(_, par_sub, on = :date, makeunique = true)
        @transform begin
            :temperature_sum = yearly_temp_cumsum(clim_sub) .* u"°C"
            :temperature = :temperature .* u"°C"
            :soiltemperature = :soiltemperature .* u"°C"
            :precipitation = :precipitation .* u"mm / d"
            :PET = :PET .* u"mm / d"
            :PAR = :PAR .* 10000 .* u"MJ / (d * ha)"
            :mowing = prepare_mowing(mow_sub) ./ 100 .* u"m"
            :grazing = prepare_grazing(graz_sub) ./ u"ha"
        end
        @select(:date,
            :temperature, :temperature_sum, :soiltemperature,
            :precipitation, :PET,
            :PAR,
            :mowing, :grazing)
    end

    daily_data = (; zip(Symbol.(names(daily_data_prep)), eachcol(daily_data_prep))...)

    ### ----------------- initial biomass
    initbiomass = 1500u"kg / ha"
    if !constant_init_biomass
        initbiomass_sub = @subset data.input.initbiomass :plotID.==plotID
        initbiomass = initbiomass_sub.biomass_init[1] * u"kg / ha"
    end

    ### ----------------- abiotic
    nut_sub = @subset data.input.nut :plotID.==plotID
    nutrient_index = nut_sub.n_index[1]
    total_N = nut_sub.Total_N[1]
    CN_ratio = nut_sub.CN_ratio[1]

    soil_sub = @subset data.input.soil :plotID.==plotID
    WHC = soil_sub.WHC[1] * u"mm"
    PWP = soil_sub.PWP[1] * u"mm"
    Clay = soil_sub.Clay[1]
    Silt = soil_sub.Silt[1]
    Sand = soil_sub.Sand[1]
    organic = soil_sub.organic[1]
    bulk = soil_sub.bulk[1]
    root_depth = soil_sub.root_depth[1]

    #### -------------- whether parts of the simulation are included
    included = (;
        senescence_included,
        potgrowth_included,
        mowing_included,
        grazing_included,
        below_included,
        height_included,
        water_red,
        nutrient_red,
        temperature_red,
        season_red,
        radiation_red)

    return (nspecies,
        npatches,
        patch_xdim = Int(sqrt(npatches)),
        patch_ydim = Int(sqrt(npatches)),
        nutheterog,
        included,
        startyear,
        endyear,
        site = (;
            initbiomass,
            nutrient_index = fill(nutrient_index, npatches),
            total_N,
            CN_ratio,
            WHC = fill(WHC, npatches),
            PWP = fill(PWP, npatches),
            Clay,
            Silt,
            Sand,
            organic,
            bulk,
            root_depth),
        daily_data,
        constant_seed)
end

function yearly_temp_cumsum(d::DataFrame)
    unqiue_years = unique(d.year)
    yearly_cumsum = Float64[]

    for y in unqiue_years
        d_sub = @subset d d.year.==y
        append!(yearly_cumsum, cumsum(d_sub.temperature))
    end

    return yearly_cumsum
end

function prepare_mowing(d::DataFrame)
    startyear = minimum(d.year)
    endyear = maximum(d.year)

    days = Dates.Date(startyear):Dates.lastdayofyear(Dates.Date(endyear))
    t = TimeArray((date = days,
            index = 1:length(days)),
        timestamp = :date)

    cutheights = d.CutHeight_cm1[.!ismissing.(d.CutHeight_cm1)]
    mean_cutheight = 7.0
    if !isempty(cutheights)
        mean_cutheight = mean(cutheights)
    end

    cutHeight_vec = zeros(length(days))

    for row in eachrow(d)
        y = row.year
        for i in 1:5
            if !ismissing.(row["MowingDay$i"])
                mowing_date = Dates.Date(y) + Dates.Day(row["MowingDay$i"])
                cH = row["CutHeight_cm$i"]
                if ismissing(cH)
                    cH = mean_cutheight
                end

                index = values(t[mowing_date].index)[1]
                cutHeight_vec[index] = cH
            end
        end
    end

    return cutHeight_vec
end

function prepare_grazing(d::DataFrame)
    startyear = minimum(d.year)
    endyear = maximum(d.year)

    days = Dates.Date(startyear):Dates.lastdayofyear(Dates.Date(endyear))
    t = TimeArray((date = days,
            index = 1:length(days)),
        timestamp = :date)

    grazing_vec = zeros(length(days))

    for row in eachrow(d)
        y = row.year
        for i in 1:4
            if !ismissing.(row["start_graz$i"])
                grazing_startdate = row["start_graz$i"]
                grazing_enddate = row["end_graz$i"]
                grazing_intensity = row["inten_graz$i"]

                if isnan(grazing_intensity)
                    @warn "Grazing intensity NaN"
                    @show grazing_startdate, grazing_enddate, d.plotID[1]
                end
                start_index = values(t[grazing_startdate].index)[1]
                end_index = values(t[grazing_enddate].index)[1]
                grazing_vec[start_index:end_index] .= grazing_intensity
            end
        end
    end
    return grazing_vec
end
