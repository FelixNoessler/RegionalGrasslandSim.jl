"""
    one_day!(;
        du_biomass,
        du_water,
        u_water,
        u_biomass,
        p,
        t)

TBW
"""
function one_day!(; calc, p, t)
    day_of_year = t - ((t - 1) ÷ 365 * 365)
    year = 1 + (t - 1) ÷ 365

    for pa in 1:(p.npatches)
        # --------------------- biomass dynamics
        ### this line is needed because it is possible
        ## that there are numerical errors
        calc.u_biomass[pa, calc.u_biomass[pa, :] .< 1e-10u"kg / ha"] .= 0u"kg / ha"
        patch_biomass = @view calc.u_biomass[pa, :]
        patch_water = calc.u_water[pa]
        calc.defoliation .= 0.0u"kg / (ha * d)"

        if any(isnan.(patch_biomass))
            @error "Patch biomass isnan: $patch_biomass" maxlog=10
        end

        if sum(patch_biomass) > 0u"kg / ha"
            # -------------- mowing
            if day_of_year ∈ p.mowing_days[year]
                mowing_index = day_of_year .== p.mowing_days[year]
                mowing_height = p.mowing_heights[year][mowing_index][1]

                last_day_index = findfirst(mowing_index) - 1
                days_since_last_mowing = 200
                if last_day_index > 0
                    days_since_last_mowing = day_of_year -
                                             p.mowing_days[year][last_day_index][1]
                end

                calc.defoliation += Growth.mowing(;
                    mowing_height,
                    days_since_last_mowing,
                    CH = p.species.CH,
                    biomass = patch_biomass,
                    mowing_mid_days = p.inf_p.mowing_mid_days)
            end

            # -------------- grazing
            calc.defoliation += Growth.grazing(;
                LD = p.grazing[t],
                biomass = patch_biomass,
                ρ = p.species.ρ,
                nspecies = p.nspecies,
                grazing_half_factor = p.inf_p.grazing_half_factor)

            calc.defoliation += Growth.trampling(;
                LD = p.grazing[t],
                biomass = patch_biomass,
                CH = p.species.CH,
                nspecies = p.nspecies,
                trampling_factor = p.inf_p.trampling_factor)

            # -------------- growth
            Growth.growth(;
                t, p, calc,
                biomass=patch_biomass,
                WR=patch_water)

            # -------------- senescence
            calc.sen = Growth.senescence(;
                ST = p.env_data.temperature_sum[t],
                biomass = patch_biomass,
                μ = p.species.μ)

        else
            @warn "Sum of patch biomass = 0" maxlog=10
        end

        # -------------- stochasticity
        # max_loss = max.(patch_biomass .+ determin .* u"d", 0u"kg / ha")
        # ds = [truncated(Normal(0, 0.5); lower=-0.2*l, upper=l) for l in ustrip.(max_loss)]
        # stochas = rand.(ds)u"kg / (ha * d)"

        # -------------- net growth
        calc.du_biomass[pa, :] = calc.act_growth - calc.sen - calc.defoliation

        # --------------------- water dynamics
        calc.du_water[pa], calc.o_evaporation[t, pa] = Water.change_water_reserve(;
            WR = calc.u_water[pa],
            precipitation = p.env_data.precipitation[t],
            LAItot = calc.LAItot,
            PET = p.env_data.PET[t],
            WHC = p.site.WHC,
            PWP = p.site.PWP)
    end

    return nothing
end

"""
    initial_conditions(p)

TBW
"""
function initial_conditions(; initbiomass, npatches, nspecies)
    u_biomass = fill(initbiomass / nspecies,
        npatches,
        nspecies)
    u_water = fill(180.0, npatches)u"mm"

    return u_biomass, u_water
end

"""
    initialize_parameters(; input_obj)

TBW
"""
function initialize_parameters(; input_obj)
    #--------- calculate leaf senescence rate μ
    sla = ustrip.(uconvert.(u"cm^2 / g", input_obj.traits.SLA))
    LL = 10 .^ ((log10.(sla) .- 2.41) ./ -0.38) .* 365.25 ./ 12 .* u"d"
    sen_intercept = input_obj.inf_p.senescence_intercept / 1000
    sen_rate = input_obj.inf_p.senescence_rate / 1000
    μ = sen_intercept * u"d^-1" .+ sen_rate .* 1 ./ LL

    #--------- grazing parameter ρ
    ρ = ustrip.(input_obj.traits.LNCM)

    #--------- functional response
    myco_nutr_right_bound, myco_nutr_midpoint = FunctionalResponse.amc_nut_response(;
        mycorrhizal_colon = input_obj.traits.AMC)

    srsa_right_bound, srsa_midpoint = FunctionalResponse.srsa_response(;
        SRSA_above = ustrip.(input_obj.traits.SRSA_above))

    sla_water_midpoint = FunctionalResponse.sla_water_response(;
        SLA = ustrip.(input_obj.traits.SLA))

    #--------- Distance matrix for below ground competition
    rel_t = input_obj.relative_traits
    trait_similarity = Growth.similarity_matrix(;
        scaled_traits = hcat(rel_t.SLA, rel_t.SRSA_above, rel_t.AMC))

    #--------- store everything in one object
    p = (;
        npatches = input_obj.npatches,
        nspecies = input_obj.nspecies,
        species = (;
            SLA = input_obj.traits.SLA,
            μ,
            ρ,
            LL,
            AMC = input_obj.traits.AMC,
            SRSA_above = input_obj.traits.SRSA_above,
            LA = input_obj.traits.LA,
            CH = input_obj.traits.CH,
            LNCM = input_obj.traits.LNCM,
            fun_response = (;
                myco_nutr_right_bound,
                myco_nutr_midpoint,
                srsa_right_bound,
                srsa_midpoint,
                sla_water_midpoint)),
        inf_p = input_obj.inf_p,
        site = input_obj.site,
        trait_similarity,
        env_data = input_obj.env_data,
        mowing_days = input_obj.mowing_days,
        mowing_heights = input_obj.mowing_heights,
        grazing = input_obj.grazing,
        water_reduction = input_obj.water_reduction,
        nutrient_reduction = input_obj.nutrient_reduction)

    return p
end


"""
    solve_prob(; input_obj)

TBW
"""
function solve_prob(; input_obj)
    p = initialize_parameters(; input_obj)
    tmax = input_obj.nyears * 365
    ts = collect(1:tmax)

    #### initial conditions of the state variables
    u_biomass, u_water = initial_conditions(;
        nspecies = p.nspecies,
        npatches = p.npatches,
        initbiomass = p.site.initbiomass)

    calc = ComponentVector(;
        ### vectors that store the state variables
        u_biomass,
        u_water,

        #### vectors that store the change of the state variables
        du_biomass = zeros(p.npatches, p.nspecies)u"kg / (ha * d)",
        du_water = zeros(p.npatches)u"mm / d",

        #### vectors that are used in the calculations
        defoliation = fill(NaN, p.nspecies)u"kg / (ha * d)",
        act_growth = fill(NaN, p.nspecies)u"kg / (ha * d)",
        pot_growth = fill(NaN, p.nspecies)u"kg / (ha * d)",
        sen = fill(NaN, p.nspecies)u"kg / (ha * d)",
        LAIs = fill(NaN, p.nspecies),
        LAItot = NaN,

        CHinfluence = fill(NaN, p.nspecies),
        Waterred = fill(NaN, p.nspecies),
        Nutred = fill(NaN, p.nspecies),
        below = fill(NaN, p.nspecies),
        species_specific_red = fill(NaN, p.nspecies),

        #### output vectors of the state variables
        o_biomass = fill(NaN, length(ts), p.npatches, p.nspecies)u"kg/ha",
        o_water = fill(NaN, length(ts), p.npatches)u"mm",

        #### outputs that are not states
        o_evaporation = fill(NaN, length(ts), p.npatches)u"mm/d"
    )

    main_loop!(; calc, ts, p)
    calc.o_biomass[iszero.(calc.o_biomass)] .= 0u"kg / ha"

    return (;
        biomass = calc.o_biomass,
        water = calc.o_water,
        evaporation = calc.o_evaporation,
        t = ts,
        p = p)
end


function main_loop!(; calc, ts, p)
    for t in ts
        one_day!(; calc, p, t)

        calc.u_biomass += calc.du_biomass * u"d"
        calc.u_water += calc.du_water * u"d"

        calc.o_biomass[t, :, :] = calc.u_biomass
        calc.o_water[t, :] = calc.u_water
    end
end
