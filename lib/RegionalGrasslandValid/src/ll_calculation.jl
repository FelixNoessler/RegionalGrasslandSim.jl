function loglikelihood_model(sim::Module; inf_p, plotID)
    inf_p = (; inf_p...)
    data, sol = get_plottingdata(sim;
        inf_p,
        plotID,
        startyear = 2012,
        endyear = 2021)

    ########################## Calculate likelihood
    ################## satellite biomass
    ### select the days where we have biomass estimated by the satellite images
    sim_biomass = @view sol.biomass[data.biomass_t, :, :]

    ### calculate the sum of biomass of all plant species
    biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

    ### calculate the likelihood
    biomass_d = MvNormal(biomass_sum, inf_p.sigma_biomass * I)
    ll_satbiomass = logpdf(biomass_d, data.biomass)

    ################## measured biomass
    f = 5 * 365 .< data.measured_biomass_t .< 10 * 365
    sim_biomass = @view sol.biomass[data.measured_biomass_t[f], :, :]

    ### calculate the sum of biomass of all plant species
    biomass_sum = vec(sum(ustrip.(sim_biomass); dims = 3))

    ### calculate the likelihood
    biomass_d = MvNormal(biomass_sum, inf_p.sigma_biomass * I)
    ll_measuredbiomass = logpdf(biomass_d, data.measured_biomass[f])

    ################## soil moisture
    sim_soilwater = ustrip.(sol.water)[data.soilmoisture_t]
    sim_soilmoisture = sim_soilwater ./ sol.p.site.root_depth .* inf_p.moisture_conv
    soilmoisture_d = MvNormal(sim_soilmoisture, inf_p.sigma_soilmoisture * I)
    ll_soilmoisture = logpdf(soilmoisture_d, data.soilmoisture)

    ################## total log likelihood
    @info """
          sat: $ll_satbiomass, mea: $ll_measuredbiomass, mois: $ll_soilmoisture
          """ maxlog=30
    ll = ll_satbiomass + ll_measuredbiomass + ll_soilmoisture

    return ll
end

VIP_plots = ["$(explo)0$i" for i in 1:9 for explo in ["HEG", "SEG", "AEG"]];

function ll_VIPS_t(sim; inf_p)
    ll = Threads.Atomic{Float64}(0.0)
    Threads.@threads for plotID in VIP_plots
        ll_plot = loglikelihood_model(sim;
            plotID,
            inf_p)

        Threads.atomic_add!(ll, ll_plot)
    end

    ### free RAM space
    GC.gc()
    ccall(:malloc_trim, Cvoid, (Cint,), 0)

    return ll[]
end

function ll_VIPS(sim; inf_p)
    ll = 0.0
    for plotID in VIP_plots
        ll_plot = loglikelihood_model(sim;
            plotID,
            inf_p)

        ll += ll_plot
    end

    return ll
end
