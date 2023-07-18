function dashboard(sim, valid, scen; datapath)
    fig = Figure(; resolution=(1700, 900))

    top_menu = fig[1,1:2] = GridLayout()
    plots_layout = fig[2,1] = GridLayout()
    param_layout = fig[2,2] = GridLayout()

    colsize!(fig.layout, 2, Relative(0.3))
    rowgap!(fig.layout, 1, 0)

    run_button = Button(top_menu[1,1]; label="run")
    left_layout = top_menu[1,2] = GridLayout()
    mowing_layout = top_menu[1,3] = GridLayout()
    center_layout = top_menu[1,4] = GridLayout()
    right_layout = top_menu[1,5] = GridLayout()
    [colgap!(top_menu, i, s)  for (i,s) in enumerate([80,50,50,80])]

    left_box = top_menu[1, 2:4]
    Box(left_box[1,1], alignmode=Outside(-8);
        color=(:black, 0.0))
    Box(right_layout[1:2,1:2], alignmode=Outside(-8);
        color=(:black, 0.0))

    ############# Exploratory
    Label(left_layout[1,1:2], "Run Scenario: Exploratory";
        tellwidth=false, halign=:left)
    menu_explo = Menu(left_layout[2,1:2], options = zip([
        "Schorfheide-Chorin",
        "Hainich",
        "Schwäbische Alb"
    ], ["SCH", "HAI", "ALB"]),
        tellwidth=false
    )


    ############# Number of years
    Label(left_layout[3,1], "Number of years";
        tellwidth=false, halign=:right)
    nyears = Observable(10)
    tb_years = Textbox(left_layout[3, 2],
        placeholder=string(nyears.val),
        stored_string=string(nyears.val),
        validator=Int64)

    on(tb_years.stored_string) do s
        nyears[] = parse(Int64, s)
    end

    ############# Number of species
    Label(left_layout[4,1], "Number of species";
        tellwidth=false, halign=:right)
    nspecies = Observable(25)
    tb_species = Textbox(left_layout[4, 2],
        placeholder=string(nspecies.val),
        stored_string=string(nspecies.val),
        validator=Int64)

    on(tb_species.stored_string) do s
        nspecies[] = parse(Int64, s)
    end

    ###########
    Label(left_layout[5,1], "include water reduction?";
        tellwidth=true, halign=:right)
    toggle_water_red = Toggle(
        left_layout[5, 2],
        active = true
    )

    Label(left_layout[6,1], "include nutrient reduction?";
        tellwidth=true, halign=:right)
    toggle_nutr_red = Toggle(
        left_layout[6, 2],
        active = true
    )

    [rowgap!(left_layout, i, dist) for (i, dist) in enumerate([5,15,4,5,5])]


    ############# Mowing
    Label(mowing_layout[1,1:2], "Mowing (date)";
        tellwidth=true, halign=:left)
    toggles_mowing = [
        Toggle(mowing_layout[i+1, 1],
            active = false, halign=:left) for i in 1:5]
    tb_mowing_date = [Textbox(mowing_layout[i+1, 2],
        halign=:left,
        placeholder=mow_date, validator=test_date,
        stored_string=mow_date) for (i,mow_date) in enumerate(
            ["05-01","06-01", "07-01", "08-01", "09-01"])
    ]
    [rowgap!(mowing_layout, i, 4) for i in 1:5]
    colgap!(mowing_layout, 1, 10)

    ############# Grazing
    Label(center_layout[1,1:4], "Grazing (start date, end date, LD)";
        halign=:left)
    toggles_grazing = [
        Toggle(center_layout[i+1, 1], active = false) for i in 1:2]
    tb_grazing_start = [Textbox(center_layout[i+1, 2],
        placeholder=graz_start, validator=test_date,
        stored_string=graz_start) for (i,graz_start) in enumerate(["05-01","08-01"])]
    tb_grazing_end = [Textbox(center_layout[i+1, 3],
        placeholder=graz_end, validator=test_date,
        stored_string=graz_end) for (i,graz_end) in enumerate(["07-01","10-01"])]
    tb_grazing_intensity = [Textbox(center_layout[i+1, 4],
        placeholder=string(intensity), validator=Float64,
        stored_string=string(intensity)) for (i,intensity) in enumerate([2, 2])]

    ############# Nutrients
    slidergrid_nut = SliderGrid(
        center_layout[4,1:4],
        (label = "Nutrient\nindex", range = 0.0:0.01:1.0, format = "{:.1f}", startvalue = 0.8))
    slider_nut = slidergrid_nut.sliders[1]

    ############# water holding capacity and permanent wilting point
    Label(center_layout[5,1], "PWP, WHC";
        halign=:left)
    slider_pwp_whc = IntervalSlider(center_layout[5, 2:3], range = LinRange(50, 250, 500),
        startvalues = (80, 150))
    whc_pwp_labeltext = lift(slider_pwp_whc.interval) do pwp_whc
        string(Int.(round.(pwp_whc)))
    end
    Label(center_layout[5, 4], whc_pwp_labeltext)
    [colgap!(center_layout, i, dist) for (i, dist) in enumerate([5,4,4])]
    [rowgap!(center_layout, i, dist) for (i, dist) in enumerate([5,4,10,10])]

    ############# Plot ID
    Label(right_layout[1,1], "OR: Plot ID";
        tellwidth=true, halign=:left)
    toggle_plotID = Toggle(
        right_layout[1, 2],
        active = false,
        tellwidth=false,
        halign=:left
    )
    menu_plotID = Menu(right_layout[2,1:2];
        options=["$(explo)0$i" for i in 1:9 for explo in ["HEG", "SEG", "AEG"]])

    ############# Abiotic variable
    Label(right_layout[3,1:2], "Abiotic variable (right lower plot)";
        tellwidth=false, halign=:left)
    menu_abiotic = Menu(right_layout[4,1:2], options = zip([
            "Precipitation [mm d⁻¹]",
            "Potential evapo-\ntranspiration [mm d⁻¹]",
            "Air temperature [°C]\n",
            "Air temperaturesum [°C]\n",
            "Photosynthetically active\nradiation [MJ m⁻² d⁻¹]"
        ], [
            :precipitation,
            :PET,
            :temperature,
            :temperature_sum,
            :PAR
            ])
        )

    ############# Coloring -> traits or biomass
    Label(right_layout[5,1:2], "Color/trait (right upper plot)";
        tellwidth=false, halign=:left)
    menu_color = Menu(right_layout[6,1:2], options = zip([
        "Specific leaf area [m² g⁻¹]",
        "Leaf nitrogen per leaf mass [mg g⁻¹]",
        "Height [m]",
        "Leaf life span [d]",
        "Mycorrhizal colonisation",
        "Root surface area /\nabove ground biomass [m² g⁻¹]",
        # "Mowing effect λ\n(part that is removed)", :λ,
        "Biomass"
    ], [:SLA, :LNCM, :CH, :LL, :AMC, :SRSA_above, :biomass])
    )
    [rowgap!(right_layout, i, dist) for (i, dist) in enumerate([5,20,5,15,5])]

    ############# Parameter values
    parameter_names = [
        "sigma_biomass",
        "sigma_evaporation",
        "sigma_soilmoisture",
        "moisture_conv",
        "senescence_intercept",
        "senescence_rate",
        "below_competition_strength",
        "trampling_factor",
        "grazing_half_factor",
        "mowing_mid_days",
        "max_SRSA_water_reduction",
        "max_SLA_water_reduction",
        "max_AMC_nut_reduction",
        "max_SRSA_nut_reduction"
    ]
    p_lower = [0.1,    0.1, 0.1, 0.1, 1e-6,  0,  0.0,   50.0,  500,  10, 0.0, 0.0, 0.0, 0.0]
    p_upper = [100000, 100, 100, 1.2, 1e-1,  10, 0.01,  200.0, 5000,  150, 1.0, 1.0, 1.0, 1.0]
    p_start = [
        82649.94, 9.590659, 0.0144827, 0.5530409, 0.07826682, 3.755483, 0.00347689,
        184.2303, 2911.013, 64.26925, 0.9826136, 0.919179, 0.812613, 0.3877618, 3446.606, 3613.506, -166.9]
    param_slider_prep = [(
        label=parameter_names[i],
        range=p_lower[i]:0.0001:p_upper[i],
        format = "{:.4f}",
        startvalue=p_start[i] ) for i in 5:length(parameter_names)
    ]
    sliders_param = SliderGrid(
        param_layout[1, 1],
        tellheight=false,
        param_slider_prep...;)

    #############
    axes = [
        Axis(
            plots_layout[i,u];
            # backgroundcolor=(:grey, i+(u-1)*3 .∈ Ref([1,4]) ? 0.1 : 0.0),
            alignmode = Outside()
            )
            for (i,u) in zip([1,1,2,2,2], [1,2,1,2,3])
    ]

    cb = Colorbar(plots_layout[1,3];
        colorrange=(-0.5, 0.5),
        valign=:bottom,
        alignmode = Inside(),
        vertical = false,
        flipaxis = true,
        tellwidth=false,
        tellheight=false,
        height=20,
        width=300)

    ###########
    sol = nothing
    plotID = nothing

    still_running = false


    on(run_button.clicks) do n

        if !still_running
            still_running = true

            # ------------- parameter values
            parameter_vals = [s.value.val for s in sliders_param.sliders]
            inf_p = (; zip(Symbol.(parameter_names[5:end]), parameter_vals)...)

            use_simulated_data = !toggle_plotID.active.val
            input_obj = nothing

            if use_simulated_data
                # ------------- mowing
                mowing_selected = [toggle.active.val for toggle in toggles_mowing]
                mowing_dates = [tb.stored_string.val for tb in tb_mowing_date][mowing_selected]

                ## final vectors of vectors
                day_of_year_mowing = Dates.dayofyear.(Dates.Date.(mowing_dates, "mm-dd"))
                mowing_days = [day_of_year_mowing for _ in 1:nyears.val]
                mowing_heights = deepcopy(mowing_days)
                [mowing_heights[n][:] .= 7 for n in 1:nyears.val]

                # ------------- grazing
                grazing_selected = [toggle.active.val for toggle in toggles_grazing]
                grazing_start = [
                    tb.stored_string.val for tb in tb_grazing_start
                ][grazing_selected]
                grazing_end = [
                    tb.stored_string.val for tb in tb_grazing_end
                ][grazing_selected]
                grazing_intensity = [
                    parse(Float64, tb.stored_string.val)
                    for tb in tb_grazing_intensity
                ][grazing_selected]

                # ------------- soil nutrients
                nutrient_index = slider_nut.value.val

                # ------------- soil water properties
                PWP, WHC = slider_pwp_whc.interval.val

                water_reduction = toggle_water_red.active.val
                nutrient_reduction = toggle_nutr_red.active.val


                input_obj = scen.scenario_input(;
                    datapath,
                    inf_p,
                    nyears=nyears.val,
                    nspecies=nspecies.val,
                    explo=menu_explo.selection.val,
                    mowing_heights,
                    mowing_days,
                    nutrient_index,
                    PWP,
                    WHC,
                    grazing_start,
                    grazing_end,
                    grazing_intensity,
                    water_reduction,
                    nutrient_reduction)
            else
                plotID = menu_plotID.selection.val
                input_obj = valid.validation_input(;
                    plotID,
                    nyears=10,
                    inf_p
                )
            end

            sol = sim.solve_prob(; input_obj)

            update_plots(; sol, menu_color, menu_abiotic, axes, cb, plotID, valid)
            still_running = false
        end
    end

    on(menu_color.selection) do n
        if !isnothing(sol)
            update_plots(; sol, menu_color, menu_abiotic, axes, cb, plotID, valid)
        end
    end

    on(menu_abiotic.selection) do n
        if !isnothing(sol)
            update_plots(; sol, menu_color, menu_abiotic, axes, cb, plotID, valid)
        end
    end



    run_button.clicks[] = 1

    return fig
end


function update_plots(; sol, menu_color, menu_abiotic, axes, cb, plotID, valid)
    trait = menu_color.selection.val
    name_index = getindex.([menu_color.options.val...], 2) .== trait
    trait_name = first.([menu_color.options.val...])[name_index][1]

    color = nothing
    if trait == :biomass
        color = ustrip.(species_biomass(sol.biomass))
    else
        color = ustrip.(sol.p.species[trait])
    end

    colormap = :viridis
    colorrange = (minimum(color), maximum(color))

    cb.colormap = colormap
    cb.colorrange = colorrange
    cb.label = trait_name

    ########### Biomass
    empty!(axes[1])
    band_patch(axes[1]; patch=1, sol, color, colormap, colorrange)

    if ! isnothing(plotID)
        data = valid.get_validation_data(; plotID)
        inityears = 5
        mdata_t = (inityears * 365 .+ data.measured_biomass_t) ./ 365
        f = mdata_t .<= 10
        mdata = data.measured_biomass ./ 1000
        mdata1 = data.measured_biomass1 ./ 1000

        scatter!(axes[1], mdata_t[f], mdata[f], marker=:hexagon, color=:black, markersize=15)
        scatter!(axes[1], mdata_t[f], mdata1[f], marker=:hexagon, color=:black, markersize=15)

    end

    ########### Growth rates
    empty!(axes[4])
    growth_rates(axes[4]; patch=1, sol, color, colormap, colorrange)


    ########### Main biomass vs traits
    empty!(axes[2])
    if menu_color.selection.val != :biomass
        trait_mean_biomass(trait, trait_name, axes[2]; sol, color, colormap, colorrange)
    else
        axes[2].xlabel = ""
        autolimits!(axes[2])
    end

    ########### Soil water
    empty!(axes[3])
    my_set = Dict(:markersize=>4, :linewidth=>0.1)
    scatterlines!(axes[3], sol.t ./ 365, ustrip.(sol.water[:, 1]);
        color=:turquoise3,
        my_set...)
    PWP = ustrip(sol.p.site.PWP)
    WHC = ustrip(sol.p.site.WHC)
    lines!(axes[3], [sol.t[1] / 365, sol.t[end] / 365], [PWP, PWP];
        color=:blue)
    lines!(axes[3], [sol.t[1] / 365, sol.t[end] / 365], [WHC, WHC];
        color=:blue)
    axes[3].ylabel="Soil water [mm]"
    axes[3].xlabel="Time [years]"
    ylims!(axes[3], 0.0, nothing)


    ########### Abiotic plot
    abiotic_colors = [:blue, :brown, :red, :red, :orange]
    abiotic = menu_abiotic.selection.val
    name_index = getindex.([menu_abiotic.options.val...], 2) .== abiotic
    abiotic_name = first.([menu_abiotic.options.val...])[name_index][1]
    abiotic_color = abiotic_colors[name_index][1]

    empty!(axes[5])
    scatterlines!(axes[5], sol.t[2:end] ./ 365, ustrip.(sol.p.env_data[abiotic]);
        color=abiotic_color,
        my_set...)
    axes[5].ylabel=abiotic_name
    axes[5].xlabel="Time [years]"

    return nothing
end
