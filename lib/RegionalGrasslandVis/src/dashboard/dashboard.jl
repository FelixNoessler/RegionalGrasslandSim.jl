function dashboard(; sim::Module, valid::Module, scen::Module, inf_p_start)
    plot_obj = dashboard_layout(; inf_p_start)

    still_running = false
    sol = nothing
    valid_data = nothing
    updating_slider = false

    on(plot_obj.obs.run_button.clicks) do n
        if !still_running
            still_running = true

            input_obj = prepare_input(; plot_obj, inf_p_start, valid, scen)
            sol = sim.solve_prob(; input_obj)
            valid_data = get_valid_data(;
                plot_obj, startyear = Dates.year(sol.date[1]), valid)
            update_plots(; sol, plot_obj, valid_data)

            updating_slider = true
            plot_obj.obs.slider_time.range = sol.t[1]:sol.t[end]
            set_close_to!(plot_obj.obs.slider_time, sol.t[end])
            updating_slider = false

            ll_obj = valid.loglikelihood_model(sim;
                inf_p = input_obj.inf_p,
                plotID = plot_obj.obs.menu_plotID.selection.val,
                nspecies = input_obj.nspecies,
                data = valid_data,
                sol = sol,
                return_seperate = true)

            plot_obj.obs.lls.biomass[] = ll_obj.biomass
            plot_obj.obs.lls.traits[] = ll_obj.trait
            plot_obj.obs.lls.soilmoisture[] = ll_obj.soilmoisture

            still_running = false
        end
    end

    on(plot_obj.obs.menu_plotID.selection) do n
        plot_obj.obs.run_button.clicks[] = 1
    end

    plot_obj.obs.run_button.clicks[] = 1

    on(plot_obj.obs.toggle_grazmow.active) do n
        band_patch(; plot_obj, sol, valid_data)
    end

    on(plot_obj.obs.toggle_bands.active) do n
        band_patch(; plot_obj, sol, valid_data)
    end

    on(plot_obj.obs.menu_color.selection) do n
        band_patch(; plot_obj, sol, valid_data)
        trait_time_plot(; plot_obj, sol, valid_data)
        trait_mean_biomass(; sol, plot_obj, t = sol.t[end])
    end

    on(plot_obj.obs.menu_abiotic.selection) do n
        abiotic_plot(; sol, plot_obj)
    end

    on(plot_obj.obs.toggle_traitvar.active) do n
        trait_time_plot(; plot_obj, sol, valid_data)
    end

    on(plot_obj.obs.toggle_validdata.active) do n
        valid_data = get_valid_data(;
            plot_obj, startyear = Dates.year(sol.date[1]), valid)
        band_patch(; plot_obj, sol, valid_data)
        trait_time_plot(; plot_obj, sol, valid_data)
        soilwater_plot(; sol, valid_data, plot_obj)
    end

    time_slider_changed = Dates.now()
    on(plot_obj.obs.slider_time.value) do t
        ms_since_change = Dates.value(Dates.now() - time_slider_changed)
        if !updating_slider && ms_since_change > 100
            updating_slider = true
            time_slider_changed = Dates.now()
            t = plot_obj.obs.slider_time.value.val
            trait_mean_biomass(; sol, plot_obj, t = t)
            updating_slider = false
        end
    end
end

function get_valid_data(; plot_obj, startyear, valid)
    is_validation = plot_obj.obs.toggle_plotID.active.val
    plotID = plot_obj.obs.menu_plotID.selection.val
    show_validdata = plot_obj.obs.toggle_validdata.active.val

    data = nothing
    if is_validation .&& show_validdata
        data = valid.get_validation_data(;
            plotID, startyear)
    end
    return data
end

function update_plots(; sol, plot_obj, valid_data)
    ########### Biomass
    band_patch(;
        plot_obj,
        sol,
        valid_data)

    ########### Trait changes over time
    trait_time_plot(;
        plot_obj, sol, valid_data)

    ########### biomass vs traits
    trait_mean_biomass(;
        sol,
        plot_obj,
        t = sol.t[end])

    ########### Soil water
    soilwater_plot(; sol, valid_data, plot_obj)

    # ########### Growth rates
    # growth_rates(axes[5]; patch = 1, sol, color, colormap, colorrange, plotID)

    ########### Abiotic plot
    abiotic_plot(; sol, plot_obj)

    return nothing
end
