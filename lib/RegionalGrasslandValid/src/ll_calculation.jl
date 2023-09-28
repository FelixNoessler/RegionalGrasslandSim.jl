function loglikelihood_model(sim::Module;
    inf_p,
    plotID,
    nspecies,
    pretty_print = false,
    return_seperate = false,
    only_likelihood = false,
    include_traits = true,
    include_soilmoisture = true,
    seed = rand(1:1000000),
    data = nothing,
    sol = nothing)

    ### do we need this line?
    inf_p = (; inf_p...)
    if isnothing(data) || isnothing(sol)
        data, sol = get_plottingdata(sim;
            inf_p,
            plotID,
            nspecies = Int(nspecies),
            startyear = 2009,
            endyear = 2021,
            seed)
    end

    selected_patch = 1

    ########################################################################
    ########################################################################
    ########################## prior
    prior = 0
    if !only_likelihood
        pdist_sigma_biomass = truncated(Normal(0, 1000); lower = 0)
        prior += logpdf(pdist_sigma_biomass, inf_p.b_biomass)
    end

    ########################################################################
    ########################################################################
    ########################## Calculate likelihood

    ########################################################################
    ################## measured biomass
    ########################################################################
    simbiomass = TimeArray((;
            biomass = vec(sum(ustrip.(sol.biomass[:, selected_patch, :]); dims = 2)),
            date = sol.date), timestamp = :date)

    ## select the dates on which measured biomass values
    ## are available
    simbiomass_sub = simbiomass[timestamp(data.measured_biomass)]
    simbiomass_vals = values(simbiomass_sub.biomass)

    if any(isnan.(simbiomass_vals))
        @warn "Biomass NaN"

        return -Inf
    end

    ### calculate the likelihood
    biomass_d = Product(Laplace.(simbiomass_vals, inf_p.b_biomass))
    ll_biomass = logpdf(biomass_d, values(data.measured_biomass.biomass))

    ########################################################################
    ################## soil moisture
    ########################################################################
    ll_soilmoisture = 0
    if include_soilmoisture
        #### downweight the likelihood because there are many observations
        weight = length(data.soilmoisture.t) / 13
        sim_soilwater = ustrip.(sol.water[data.soilmoisture.t, selected_patch])
        transformed_data = @. sol.p.site.root_depth *
                              (data.soilmoisture.val - inf_p.moistureconv_alpha) /
                              inf_p.moistureconv_beta
        soilmoisture_d = Product(Laplace.(sim_soilwater, inf_p.b_soilmoisture))
        ll_soilmoisture += logpdf(soilmoisture_d, transformed_data) / weight
    end

    ########################################################################
    ################## cwm trait likelihood
    ########################################################################
    ll_trait = 0
    if include_traits
        ########## prepare biomass for CWM calculation
        patch = 1
        biomass_vals = ustrip.(sol.biomass[data.traits.t, patch, :])
        total_biomass = sum(biomass_vals, dims = 2)

        ## cannot calculate cwm trait for zero biomass
        if any(iszero.(total_biomass))
            if return_seperate
                return (;
                    biomass = ll_biomass,
                    trait = -Inf,
                    soilmoisture = ll_soilmoisture)
            end

            return -Inf
        end

        relative_biomass = biomass_vals ./ total_biomass
        ntraits = length(data.traits.dim)

        for trait_name in data.traits.dim
            ### calculate CWM
            trait_vals = ustrip.(sol.p.species[trait_name])
            weighted_trait = trait_vals .* relative_biomass'
            sim_cwm_trait = vec(sum(weighted_trait; dims = 1))

            ### "measured" traits (calculated cwm from observed vegetation)
            mtraits_vals = vec(data.traits.val[:, trait_name .== data.traits.dim])

            ### calculate the likelihood
            trait_scale = Symbol(:b_, trait_name)

            trait_d = Product(Laplace.(sim_cwm_trait, inf_p[trait_scale]))
            ll = logpdf(trait_d, mtraits_vals)
            ll_trait += ll / ntraits
        end
    end

    ########################################################################
    ################## total likelihood
    ########################################################################
    ll = ll_biomass + ll_trait + ll_soilmoisture

    ########################################################################
    ################## printing
    ########################################################################
    if pretty_print
        bl, tl, = round(ll_biomass), round(ll_trait)
        sl, pl = round(ll_soilmoisture), round(prior)
        @info "biomass: $(bl) trait: $(tl) moi: $(sl) prior: $(pl)" maxlog=1000
    end

    if return_seperate
        return (; biomass = ll_biomass, trait = ll_trait,
            soilmoisture = ll_soilmoisture)
    end

    return ll + prior
end

VIP_plots = ["$(explo)0$i" for i in 1:9 for explo in ["HEG", "SEG", "AEG"]];

function ll_VIPS_t(sim; inf_p, nspecies)
    ll = Threads.Atomic{Float64}(0.0)
    Threads.@threads for plotID in VIP_plots
        ll_plot = loglikelihood_model(sim;
            plotID,
            inf_p,
            nspecies)
        Threads.atomic_add!(ll, ll_plot)
    end

    ### free RAM space
    # GC.gc()
    # ccall(:malloc_trim, Cvoid, (Cint,), 0)

    return ll[]
end

# function ll_VIPS(sim; inf_p)
#     ll = 0.0
#     for plotID in VIP_plots
#         ll_plot = loglikelihood_model(sim;
#             plotID,
#             inf_p)

#         ll += ll_plot
#     end

#     return ll
# end
