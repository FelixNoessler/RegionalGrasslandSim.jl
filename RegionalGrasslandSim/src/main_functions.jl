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
function one_day!(;
    du_biomass,
    du_water,
    u_water,
    u_biomass,
    o_evaporation,
    defoliation,
    act_growth,
    sen,
    p,
    t)

    day_of_year = t - ((t - 1) ÷ 365 * 365)
    year = 1 + (t - 1) ÷ 365

    for pa in 1:p.npatches
        # --------------------- biomass dynamics
        ### this line is needed because it is possible
        ## that there are numerical errors
        u_biomass[pa, u_biomass[pa, :] .< 1e-10u"kg / ha"] .= 0u"kg / ha"
        patch_biomass = u_biomass[pa, :]

        defoliation .= 0.0u"kg / (ha * d)"
        act_growth .= 0.0u"kg / (ha * d)"
        sen .= 0.0u"kg / (ha * d)"
        LAI_tot = 0

        if any(patch_biomass .< 0u"kg / ha")
            @error "Patch biomass below zero: $patch_biomass"
        end

        # if any(isnan.(patch_biomass))
        #     @error "Patch biomass isnan: $patch_biomass"
        # end

        if sum(patch_biomass) > 0u"kg / ha"

            # -------------- mowing
            if day_of_year ∈ p.mowing_days[year]
                mowing_index = day_of_year .== p.mowing_days[year]
                mowing_height = p.mowing_heights[year][mowing_index][1]

                last_day_index = findfirst(mowing_index) - 1
                days_since_last_mowing = 200

                if last_day_index > 0
                    days_since_last_mowing = day_of_year - p.mowing_days[year][last_day_index][1]
                end

                mown_biomass = Growth.mowing(;
                    mowing_height,
                    days_since_last_mowing,
                    CH=p.species.CH,
                    biomass=patch_biomass,
                    mowing_mid_days=p.inf_p.mowing_mid_days)
                defoliation += mown_biomass
            end

            # -------------- grazing
            grazed_biomass = Growth.grazing(;
                LD=p.grazing[t],
                biomass=patch_biomass,
                ρ=p.species.ρ,
                nspecies=p.nspecies,
                grazing_half_factor=p.inf_p.grazing_half_factor)
            defoliation += grazed_biomass

            trampled_biomass = Growth.trampling(;
                LD=p.grazing[t],
                biomass=patch_biomass,
                CH=p.species.CH,
                nspecies=p.nspecies,
                trampling_factor=p.inf_p.trampling_factor)
            defoliation += trampled_biomass

            # -------------- growth
            act_growth, LAI_tot = Growth.growth(;
                nspecies=p.nspecies,
                fun_response=p.species.fun_response,
                trait_similarity=p.trait_similarity,
                below_competition_strength=p.inf_p.below_competition_strength,
                SLA=p.species.SLA,
                CH=p.species.CH,
                biomass=patch_biomass,
                PAR=p.env_data.PAR[t],
                WR=u_water[pa],
                WHC=p.site.WHC,
                PWP=p.site.PWP,
                nutrients=p.site.nutrient_index,
                PET=p.env_data.PET[t],
                T=p.env_data.temperature[t],
                ST=p.env_data.temperature_sum[t],
                water_red=p.water_reduction,
                nutrient_red=p.nutrient_reduction,
                max_SRSA_water_reduction=p.inf_p.max_SRSA_water_reduction,
                max_SLA_water_reduction=p.inf_p.max_SLA_water_reduction,
                max_AMC_nut_reduction=p.inf_p.max_AMC_nut_reduction,
                max_SRSA_nut_reduction=p.inf_p.max_SRSA_nut_reduction)

            if any(act_growth .< 0u"kg / (ha * d)")
                @error "act_growth below zero: $act_growth"
            end
            if LAI_tot < 0
                @error "LAI_tot below zero: $LAI_tot"
            end

            # -------------- senescence
            sen = Growth.senescence(;
                                    ST=p.env_data.temperature_sum[t],
                                    biomass=patch_biomass,
                                    μ = p.species.μ)

        else
            # @info "Sum of patch biomass !> 0"
        end


        # -------------- deterministic
        determin = act_growth .- sen .- defoliation

        # -------------- stochasticity
        # max_loss = max.(patch_biomass .+ determin .* u"d", 0u"kg / ha")
        # ds = [truncated(Normal(0, 0.5); lower=-0.2*l, upper=l) for l in ustrip.(max_loss)]
        # stochas = rand.(ds)u"kg / (ha * d)"

        # -------------- net growth
        du_biomass[pa, :] = determin #.+ stochas

        # --------------------- water dynamics
        water_change, evaporation = Water.change_water_reserve(;
            WR=u_water[pa],
            precipitation=p.env_data.precipitation[t],
            LAI_tot,
            PET=p.env_data.PET[t],
            WHC=p.site.WHC,
            PWP=p.site.PWP)

        o_evaporation[pa] = evaporation
        du_water[pa] = water_change

    end

    return nothing
end

"""
    initial_conditions(p)

TBW
"""
function initial_conditions(; npatches, nspecies)
    u_biomass = fill(
        500 / nspecies,
        npatches,
        nspecies)u"kg / ha"
    u_water = fill(180.0, npatches)u"mm"

    return u_biomass, u_water
end

function similarity_matrix(; traits)
    nspecies, ntraits = size(traits)
    similarity_mat = Array{Float64}(undef, nspecies, nspecies)

    for i in 1:nspecies
        diff = zeros(nspecies)

        for n in 1:ntraits
            species_trait = traits[i, n]
            diff .+= abs.(species_trait .- traits[: , n])
        end

        similarity_mat[i, :] .= 1 .- (diff ./ ntraits)
    end

    return similarity_mat
end


"""
    initialize_parameters(; input_obj)

TBW
"""
function initialize_parameters(; input_obj)
    #--------- calculate leaf senescence rate μ
    sla = ustrip.(uconvert.(u"cm^2 / g", input_obj.traits.SLA))
    LL = 10 .^ ( (log10.(sla) .- 2.41) ./ -0.38) .* 365.25 ./ 12 .* u"d"
    sen_intercept = input_obj.inf_p.senescence_intercept/1000
    sen_rate = input_obj.inf_p.senescence_rate/1000
    μ = sen_intercept * u"d^-1" .+ sen_rate .* 1 ./ LL

    #--------- grazing parameter ρ
    ρ = ustrip.(input_obj.traits.LNCM)

    #--------- functional response
    myco_nutr_right_bound, myco_nutr_midpoint =
        FunctionalResponse.amc_nut_response(;
            mycorrhizal_colon=input_obj.traits.AMC)

    srsa_right_bound, srsa_midpoint =
        FunctionalResponse.srsa_response(;
        SRSA_above=ustrip.(input_obj.traits.SRSA_above))

    sla_water_midpoint =
        FunctionalResponse.sla_water_response(;
            SLA=ustrip.(input_obj.traits.SLA)
    )

    #--------- Distance matrix for below ground competition
    rel_t = input_obj.relative_traits
    trait_similarity = similarity_matrix(;
        traits=hcat(rel_t.SLA, rel_t.SRSA_above, rel_t.AMC))

    #--------- store everything in one object
    p = (;
         npatches=input_obj.npatches,
         nspecies=input_obj.nspecies,
         species = (;
                    SLA=input_obj.traits.SLA,
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
                        sla_water_midpoint
                    )
         ),
         inf_p=input_obj.inf_p,
         site=input_obj.site,
         trait_similarity,
         env_data=input_obj.env_data,
         mowing_days=input_obj.mowing_days,
         mowing_heights=input_obj.mowing_heights,
         grazing=input_obj.grazing,
         water_reduction=input_obj.water_reduction,
         nutrient_reduction=input_obj.nutrient_reduction)

    return p
end



"""
    solve_prob(; input_obj)

TBW
"""
function solve_prob(; input_obj)
    p = initialize_parameters(; input_obj)
    tmax = input_obj.nyears * 365
    ts = collect(0:tmax)

    #### initial conditions of the state variables
    u_biomass, u_water = initial_conditions(;
        nspecies=p.nspecies,
        npatches=p.npatches)

    #### output vectors of the state variables
    u_output_biomass = Array{Quantity{Float64}}(undef, length(ts), p.npatches, p.nspecies)
    u_output_water = Array{Quantity{Float64}}(undef, length(ts), p.npatches)

    u_output_biomass[1, :, :] = u_biomass
    u_output_water[1, :] = u_water

    #### outputs that are not states
    o_evaporation = fill(NaN, p.npatches)u"mm/d"
    o_output_evaporation = Array{Quantity{Float64}}(undef, length(ts), p.npatches)
    o_output_evaporation[1, :] = o_evaporation

    #### vectors that store the change of the state variables
    du_biomass = Array{Quantity{Float64}}(undef, p.npatches, p.nspecies)
    du_water = Array{Quantity{Float64}}(undef, p.npatches)


    #### vectors that are used in the calculations
    defoliation = zeros(p.nspecies)u"kg / (ha * d)"
    act_growth = zeros(p.nspecies)u"kg / (ha * d)"
    sen = zeros(p.nspecies)u"kg / (ha * d)"

    for (i,t) in enumerate(ts[2:end])
        one_day!(;
            du_biomass, du_water,
            u_biomass, u_water,
            o_evaporation,
            defoliation,
            act_growth,
            sen,
            p, t)

        u_biomass += du_biomass * u"d"
        u_water += du_water * u"d"

        u_output_biomass[i+1, :, :] = u_biomass
        u_output_water[i+1, :] = u_water

        o_output_evaporation[i+1, :] = o_evaporation
    end

    u_output_biomass[iszero.(u_output_biomass)] .= 0u"kg / ha"

    return (;
        biomass=u_output_biomass,
        water=u_output_water,
        evaporation=o_output_evaporation,
        t=ts,
        p=p)
end
