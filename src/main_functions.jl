"""
    one_day!(; calc, p, t)

Calculate the density differences of all state variables of one day.
"""
function one_day!(; calc, p, t)
    LAItot = 0.0

    for pa in 1:(p.npatches)
        # --------------------- biomass dynamics
        ### this line is needed because it is possible
        ## that there are numerical errors
        patch_biomass = @view calc.u_biomass[pa, :]
        calc.very_low_biomass .= patch_biomass .< 1e-30u"kg / ha" .&&
                                 .!iszero.(patch_biomass)
        if any(calc.very_low_biomass)
            patch_biomass[calc.very_low_biomass] .= 0u"kg / ha"
            @warn "Set some values of u_biomass to zero" maxlog=10
        end

        patch_water = calc.u_water[pa]
        calc.defoliation .= 0.0u"kg / (ha * d)"

        calc.nan_biomass .= isnan.(patch_biomass)
        if any(calc.nan_biomass)
            @error "Patch biomass isnan: $patch_biomass" maxlog=10
        end

        if !iszero(sum(patch_biomass))
            if p.included.mowing_included
                # -------------- mowing
                mowing_height = p.daily_data.mowing[t]
                if !iszero(mowing_height)
                    tstart = t - 200
                    tstart = tstart < 1 ? 1 : tstart
                    mowing_last200 = @view p.daily_data.mowing[tstart:t]
                    days_since_last_mowing = maximum(201 .-
                                                     findall(.!iszero.(mowing_last200)))
                    days_since_last_mowing = iszero(days_since_last_mowing) ? 200 :
                                             days_since_last_mowing

                    Growth.mowing!(;
                        calc,
                        mowing_height,
                        days_since_last_mowing,
                        height = p.species.height,
                        biomass = patch_biomass,
                        mowing_mid_days = p.inf_p.mowing_mid_days)
                end
            end

            if p.included.grazing_included
                LD = p.daily_data.grazing[t]
                if !iszero(LD)
                    # -------------- grazing
                    Growth.grazing!(;
                        calc,
                        LD,
                        biomass = patch_biomass,
                        ρ = p.species.ρ,
                        grazing_half_factor = p.inf_p.grazing_half_factor)

                    Growth.trampling!(;
                        calc,
                        LD,
                        biomass = patch_biomass,
                        height = p.species.height,
                        trampling_factor = p.inf_p.trampling_factor)
                end
            end

            # -------------- growth
            LAItot = Growth.growth!(;
                t, p, calc,
                biomass = patch_biomass,
                WR = patch_water)

            # -------------- senescence
            if p.included.senescence_included
                Growth.senescence!(;
                    sen = calc.sen,
                    ST = p.daily_data.temperature_sum[t],
                    biomass = patch_biomass,
                    μ = p.species.μ)
            end

        else
            @warn "Sum of patch biomass = 0" maxlog=10
        end

        # -------------- stochasticity
        # max_loss = max.(patch_biomass .+ determin .* u"d", 0u"kg / ha")
        # ds = [truncated(Normal(0, 0.5); lower=-0.2*l, upper=l) for l in ustrip.(max_loss)]
        # stochas = rand.(ds)u"kg / (ha * d)"

        # -------------- net growth
        calc.du_biomass[pa, :] .= calc.act_growth .- calc.sen .- calc.defoliation

        # --------------------- water dynamics
        ## no allocations!
        calc.du_water[pa] = Water.change_water_reserve(;
            WR = calc.u_water[pa],
            precipitation = p.daily_data.precipitation[t],
            LAItot,
            PET = p.daily_data.PET[t],
            WHC = p.site.WHC,
            PWP = p.site.PWP)
    end

    return nothing
end

"""
    initialize_parameters(; input_obj)

Initialize all parameters and store them in one object.
"""
function initialize_parameters(; input_obj, inf_p)
    #--------- create random traits
    seed = 120
    if !input_obj.constant_seed
        seed = rand(1:100000)
    end
    traits, rel_traits = Traits.generate_random_traits(input_obj.nspecies; seed)

    #--------- Distance matrix for below ground competition
    trait_similarity = Traits.similarity_matrix(;
        scaled_traits = hcat(rel_traits.SRSA_above, rel_traits.AMC),
        belowtrait_similarity_exponent = inf_p.belowtrait_similarity_exponent)

    #--------- calculate leaf senescence rate μ
    sla = ustrip.(uconvert.(u"cm^2 / g", traits.SLA))
    LL = 10 .^ ((log10.(sla) .- 2.41) ./ -0.38) .* 365.25 ./ 12 .* u"d"
    sen_intercept = inf_p.senescence_intercept / 1000
    sen_rate = inf_p.senescence_rate / 1000
    μ = sen_intercept * u"d^-1" .+ sen_rate ./ LL

    #--------- grazing parameter ρ [dimensionless]
    ρ = traits.LNCM .* u"g / mg"

    #--------- functional response
    myco_nut_lower, myco_nut_upper, myco_nut_midpoint = FunctionalResponse.amc_nut_response(;
        mycorrhizal_colon = traits.AMC,
        maximal_reduction = inf_p.max_AMC_nut_reduction)

    srsa_water_lower, srsa_water_upper, srsa_midpoint = FunctionalResponse.srsa_response(;
        SRSA_above = traits.SRSA_above,
        maximal_reduction = inf_p.max_SRSA_water_reduction)

    srsa_nut_lower, srsa_nut_upper, srsa_midpoint = FunctionalResponse.srsa_response(;
        SRSA_above = traits.SRSA_above,
        maximal_reduction = inf_p.max_SRSA_nut_reduction)

    sla_water_lower, sla_water_midpoint = FunctionalResponse.sla_water_response(;
        SLA = traits.SLA,
        maximal_reduction = inf_p.max_SLA_water_reduction)

    #--------- store everything in one object
    p = (;
        input_obj...,
        species = (;
            SLA = traits.SLA,
            μ,
            ρ,
            LL,
            AMC = traits.AMC,
            SRSA_above = traits.SRSA_above,
            LA = traits.LA,
            height = traits.height,
            LNCM = traits.LNCM,
            fun_response = (;
                myco_nut_lower,
                myco_nut_upper,
                myco_nut_midpoint,
                srsa_water_lower,
                srsa_water_upper,
                srsa_nut_lower,
                srsa_nut_upper,
                srsa_midpoint,
                sla_water_lower,
                sla_water_midpoint)),
        inf_p,
        trait_similarity)

    return p
end

function reset_initialconditions!(; calc, input_obj)
    calc.u_biomass .= input_obj.site.initbiomass / input_obj.nspecies
    calc.u_water .= 180.0u"mm"
end

function preallocate_vectors(; input_obj)
    n = length(input_obj.daily_data.date)

    calc = (;
        ############ vectors that store the state variables
        u_biomass = fill(input_obj.site.initbiomass / input_obj.nspecies,
            input_obj.npatches,
            input_obj.nspecies),
        u_water = fill(180.0, input_obj.npatches)u"mm",

        ############ vectors that store the change of the state variables
        du_biomass = zeros(input_obj.npatches, input_obj.nspecies)u"kg / (ha * d)",
        du_water = zeros(input_obj.npatches)u"mm / d",

        ############ output vectors of the state variables
        o_biomass = fill(NaN, n, input_obj.npatches, input_obj.nspecies)u"kg/ha",
        o_water = fill(NaN, n, input_obj.npatches)u"mm",

        ############ preallaocated vectors that are used in the calculations
        pot_growth = fill(NaN, input_obj.nspecies)u"kg / (ha * d)",
        act_growth = fill(NaN, input_obj.nspecies)u"kg / (ha * d)",
        defoliation = fill(NaN, input_obj.nspecies)u"kg / (ha * d)",
        sen = zeros(input_obj.nspecies)u"kg / (ha * d)",
        species_specific_red = fill(NaN, input_obj.nspecies),
        LAIs = fill(NaN, input_obj.nspecies),

        ## warnings, debugging, avoid numerical errors
        very_low_biomass = fill(false, input_obj.nspecies),
        nan_biomass = fill(false, input_obj.nspecies),
        neg_act_growth = fill(false, input_obj.nspecies),

        ## below ground competition
        below_split = fill(NaN, input_obj.nspecies),
        traitsimilarity_biomass = fill(NaN, input_obj.nspecies),

        ## height influence
        heightinfluence = fill(NaN, input_obj.nspecies),
        biomass_height = fill(NaN, input_obj.nspecies)u"m * kg / ha",

        ## nutrient reducer function
        nutrients_splitted = fill(NaN, input_obj.nspecies),
        Nutred = fill(NaN, input_obj.nspecies),
        amc_nut = fill(NaN, input_obj.nspecies),
        srsa_nut = fill(NaN, input_obj.nspecies),

        ## water reducer function
        water_splitted = fill(NaN, input_obj.nspecies),
        Waterred = fill(NaN, input_obj.nspecies),
        sla_water = fill(NaN, input_obj.nspecies),
        srsa_water = fill(NaN, input_obj.nspecies),

        ## mowing
        mown_height = fill(NaN, input_obj.nspecies)u"m",
        mowing_λ = fill(NaN, input_obj.nspecies),

        ## grazing
        biomass_ρ = fill(NaN, input_obj.nspecies)u"kg / ha",
        grazed_share = fill(NaN, input_obj.nspecies),

        ## trampling
        trampling_ω = fill(NaN, input_obj.nspecies)u"1/ha",
        trampled_biomass = fill(NaN, input_obj.nspecies)u"kg / ha",
        trampling_high_LD = fill(false, input_obj.nspecies))

    return calc
end

"""
    solve_prob(; input_obj, inf_p, calc = nothing)

Solve the model for one site.

Intialize the parameters, the state variables and the output vectors.
In addition some vectors are preallocated to avoid allocations in the main loop.
Then, run the main loop and store the results in the output vectors.
"""
function solve_prob(; input_obj, inf_p, calc = nothing)
    p = initialize_parameters(; input_obj, inf_p)

    if isnothing(calc)
        calc = preallocate_vectors(; input_obj)
    else
        reset_initialconditions!(; calc, input_obj)
    end

    ts = 1:length(p.daily_data.date)
    main_loop!(; calc, ts, p)

    calc.o_biomass[calc.o_biomass .< 0u"kg / ha"] .= 0u"kg / ha"

    return (;
        biomass = calc.o_biomass,
        water = calc.o_water,
        t = ts,
        date = p.daily_data.date,
        numeric_date = to_numeric.(p.daily_data.date),
        p = p)
end

"""
    main_loop!(; calc, ts, p)

Run the main loop for all days.
"""
function main_loop!(; calc, ts, p)
    for t in ts
        one_day!(; calc, p, t)

        ### no allocations! :)
        calc.u_biomass .+= calc.du_biomass .* u"d"
        calc.u_water .+= calc.du_water .* u"d"
        calc.o_biomass[t, :, :] .= calc.u_biomass
        calc.o_water[t, :] .= calc.u_water
    end
end

function to_numeric(d::Dates.Date)
    daysinyear = Dates.daysinyear(Dates.year(d))
    return Dates.year(d) + (Dates.dayofyear(d) - 1) / daysinyear
end
