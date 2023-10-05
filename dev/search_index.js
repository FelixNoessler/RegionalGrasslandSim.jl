var documenterSearchIndex = {"docs":
[{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"CurrentModule=RegionalGrasslandSim","category":"page"},{"location":"Modelling_API/Functional_response/#Functional-response","page":"Functional response","title":"Functional response","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"The growth of the plants is limited by soil water and nutrients","category":"page"},{"location":"Modelling_API/Functional_response/#water_stress","page":"Functional response","title":"Water stress","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"The species differ in the response to water stress by the different specific leaf areas and root surface areas per above ground biomass. The values of both response curves are multiplied to get the growth reduction factor.","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"It is implemented in Growth.water_reduction!.","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"","category":"page"},{"location":"Modelling_API/Functional_response/#sla","page":"Functional response","title":"Specific leaf area","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"the core of the functional response is build in FunctionalResponse.sla_water_response\nthe strength of the reduction is modified by the parameter max_SLA_water_reduction in Growth.sla_water_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SLA_water_reduction equals 1: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SLA_water_reduction equals 0.5: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"FunctionalResponse.sla_water_response\nGrowth.sla_water_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.FunctionalResponse.sla_water_response","page":"Functional response","title":"RegionalGrasslandSim.FunctionalResponse.sla_water_response","text":"sla_water_response(;\n    SLA,\n    mid_SLA = 0.025u\"m^2 / g\",\n    slope_func_parameter = 75u\"g / m^2\",\n    min_SLA_half_response = -0.8,\n    max_SLA_half_response = 0.8,\n    maximal_reduction)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.Growth.sla_water_reduction!","page":"Functional response","title":"RegionalGrasslandSim.Growth.sla_water_reduction!","text":"sla_water_reduction!(;\n    sla_water,\n    fun_response,\n    x)\n\nReduction of growth due to stronger water stress for higher specific leaf area (SLA).\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"","category":"page"},{"location":"Modelling_API/Functional_response/#srsa_water","page":"Functional response","title":"Root surface area / aboveground biomass","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"the core of the functional response is build in FunctionalResponse.srsa_response\nthe strength of the reduction is modified by the parameter max_SRSA_water_reduction in Growth.srsa_water_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SRSA_water_reduction equals 1: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SRSA_water_reduction equals 0.5: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"FunctionalResponse.srsa_response\nGrowth.srsa_water_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.FunctionalResponse.srsa_response","page":"Functional response","title":"RegionalGrasslandSim.FunctionalResponse.srsa_response","text":"srsa_response(;\n    SRSA_above,\n    mid_SRSA_above = 0.12u\"m^2 / g\",\n    slope_func_parameters = 40u\"g / m^2 \",\n    min_right_upper_bound = 0.7,\n    max_right_upper_bound = 1,\n    min_SRSA_half_response = 0.05,\n    max_SRSA_half_response = 0.6,\n    maximal_reduction)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.Growth.srsa_water_reduction!","page":"Functional response","title":"RegionalGrasslandSim.Growth.srsa_water_reduction!","text":"srsa_water_reduction!(; srsa_water, fun_response, x)\n\nReduction of growth due to stronger water stress for lower specific root surface area per above ground biomass (SRSA_above).\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/#nut_stress","page":"Functional response","title":"Nutrient stress","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"The species differ in the response to nutrient availability by different proportion of mycorrhizal colonisations  and root surface per above ground biomass. The maximum of both response curves is used for the nutrient reduction function. It is assumed that the plants needs either many fine roots per above ground biomass or have a strong symbiosis with mycorrhizal fungi. ","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"It is implemented in Growth.nutrient_reduction!.","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"","category":"page"},{"location":"Modelling_API/Functional_response/#amc","page":"Functional response","title":"Arbuscular mycorrhizal colonisation","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"the core of the functional response is build in FunctionalResponse.amc_nut_response\nthe strength of the reduction is modified by the parameter max_AMC_nut_reduction in Growth.amc_nut_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_AMC_nut_reduction equals 1: (Image: Graphical overview of the AMC functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_AMC_nut_reduction equals 0.5: (Image: Graphical overview of the AMC functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"FunctionalResponse.amc_nut_response\nGrowth.amc_nut_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.FunctionalResponse.amc_nut_response","page":"Functional response","title":"RegionalGrasslandSim.FunctionalResponse.amc_nut_response","text":"amc_nut_response(;\n    mycorrhizal_colon,\n    max_right_upper_bound = 1,\n    min_right_upper_bound = 0.7,\n    max_AMC_half_response = 0.6,\n    min_AMC_half_response = 0.05,\n    mid_AMC = 0.35,\n    slope = 10,\n    maximal_reduction)\n\nTransforms the mycorrhizal colonisation into parameters of the response curve of growth in relation to nutrient availability.\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.Growth.amc_nut_reduction!","page":"Functional response","title":"RegionalGrasslandSim.Growth.amc_nut_reduction!","text":"amc_nut_reduction!(; amc_nut, fun_response, x)\n\nReduction of growth due to stronger nutrient stress for lower arbuscular mycorrhizal colonization (AMC).\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"","category":"page"},{"location":"Modelling_API/Functional_response/#srsa_nut","page":"Functional response","title":"Root surface area / aboveground biomass","text":"","category":"section"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"the core of the functional response is build in FunctionalResponse.srsa_response\nthe strength of the reduction is modified by the parameter max_SRSA_nut_reduction in Growth.srsa_nut_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SRSA_nut_reduction equals 1: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"max_SRSA_nut_reduction equals 0.5: (Image: Graphical overview of the functional response)","category":"page"},{"location":"Modelling_API/Functional_response/","page":"Functional response","title":"Functional response","text":"Growth.srsa_nut_reduction!","category":"page"},{"location":"Modelling_API/Functional_response/#RegionalGrasslandSim.Growth.srsa_nut_reduction!","page":"Functional response","title":"RegionalGrasslandSim.Growth.srsa_nut_reduction!","text":"srsa_nut_reduction!(; srsa_nut, fun_response, x)\n\nReduction of growth due to stronger nutrient stress for lower specific root surface area per above ground biomass (SRSA_above).\n\n\n\n\n\n","category":"function"},{"location":"References/#References","page":"References","title":"References","text":"","category":"section"},{"location":"References/","page":"References","title":"References","text":"<div class=\"citation canonical\"><dl><dt>[1]</dt>\n<dd>\n<div id=\"Jouven2006\">M. Jouven, P. Carrere and R. Baumont. <i>Model predicting dynamics of biomass,  structure and digestibility of herbage in managed permanent pastures. 1. Model description</i>. <a href='https://doi.org/10.1111/j.1365-2494.2006.00515.x'>Grass and Forage Science <b>61</b>, 112–124 (2006)</a>.</div>\n</dd><dt>[2]</dt>\n<dd>\n<div id=\"Schapendonk1998\">A. Schapendonk, W. Stol, D. Kraalingen and B. Bouman. <i>LINGRA,  a sink/source model to simulate grassland productivity in Europe</i>. <a href='https://doi.org/10.1016/s1161-0301(98)00027-6'>European Journal of Agronomy <b>9</b>, 87–100 (1998)</a>.</div>\n</dd><dt>[3]</dt>\n<dd>\n<div id=\"Moulin2021\">T. Moulin, A. Perasso, P. Calanca and F. Gillet. <i>DynaGraM: A process-based model to simulate multi-species plant community dynamics in managed grasslands</i>. <a href='https://doi.org/10.1016/j.ecolmodel.2020.109345'>Ecological Modelling <b>439</b>, 109345 (2021)</a>.</div>\n</dd><dt>[4]</dt>\n<dd>\n<div id=\"Gillet2008\">F. Gillet. <i>Modelling vegetation dynamics in heterogeneous pasture-woodland landscapes</i>. <a href='https://doi.org/10.1016/j.ecolmodel.2008.05.013'>Ecological Modelling <b>217</b>, 1-18 (2008)</a>.</div>\n</dd>\n</dl></div>","category":"page"},{"location":"Modelling_API/Water_dynamics/#Water-dynamics","page":"Water dynamics","title":"Water dynamics","text":"","category":"section"},{"location":"Modelling_API/Water_dynamics/","page":"Water dynamics","title":"Water dynamics","text":"Modules = [RegionalGrasslandSim.Water]\nPublic = true\nPrivate = true","category":"page"},{"location":"Modelling_API/Water_dynamics/#RegionalGrasslandSim.Water.actual_evapotranspiration-Tuple{}","page":"Water dynamics","title":"RegionalGrasslandSim.Water.actual_evapotranspiration","text":"actual_evapotranspiration(; WR, ATr, AEv)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Water_dynamics/#RegionalGrasslandSim.Water.change_water_reserve-Tuple{}","page":"Water dynamics","title":"RegionalGrasslandSim.Water.change_water_reserve","text":"change_water_reserve(;\n    WR, precipitation, LAItot,\n    PET, WHC, PWP)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Water_dynamics/#RegionalGrasslandSim.Water.evaporation-Tuple{}","page":"Water dynamics","title":"RegionalGrasslandSim.Water.evaporation","text":"evaporation(; WR, WHC, PET, LAItot)\n\nAev(t) = WR(t) / WHC * PET(t)*[1 - min(1; LAI_tot(t)/3) ]\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Water_dynamics/#RegionalGrasslandSim.Water.transpiration-Tuple{}","page":"Water dynamics","title":"RegionalGrasslandSim.Water.transpiration","text":"transpiration(;\nWR,\nPWP, WHC,\nPET,\nLAItot)\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Water_dynamics/#RegionalGrasslandSim.Water.water_drainage-Tuple{}","page":"Water dynamics","title":"RegionalGrasslandSim.Water.water_drainage","text":"water_drainage(; WR, precipitation, WHC, AET)\n\nΔ(𝑡) = max(𝑊𝑅(𝑡) + 𝑃 (𝑡) − 𝐴𝐸𝑇 (𝑡) − 𝑊𝐻𝐶 ; 0)\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Traits/#Growth","page":"Traits","title":"Growth","text":"","category":"section"},{"location":"Modelling_API/Traits/","page":"Traits","title":"Traits","text":"CurrentModule = RegionalGrasslandSim","category":"page"},{"location":"Modelling_API/Traits/#trait_similarity","page":"Traits","title":"Calculate the trait similarity","text":"","category":"section"},{"location":"Modelling_API/Traits/","page":"Traits","title":"Traits","text":"Traits.similarity_matrix","category":"page"},{"location":"Modelling_API/Traits/#RegionalGrasslandSim.Traits.similarity_matrix","page":"Traits","title":"RegionalGrasslandSim.Traits.similarity_matrix","text":"similarity_matrix(; scaled_traits, similarity_exponent)\n\nComputes the trait similarity of all plant species.\n\nThe trait similarity between plant species i and plant species u for T traits is calculated as follows:\n\ntexttrait_similarity_iu =\n    1-fracsum_t=1^t=T\n        textscaled_trait_ti - textscaled_trait_tuT\n\nTo give each functional trait an equal influence, the trait values have been scaled by the 5 % (Q_005 t) and 95 % quantile (Q_095 t) of trait values of 100 plant species:\n\ntextscaled_trait_ti =\n    fractexttrait_ti - Q_005 t\n    Q_095 t - Q_005 t\n\nIf the rescaled trait values were below zero or above one, the values were set to zero or one respectively.\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#Growth","page":"Growth","title":"Growth","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"the net growth of the plants is modelled by...","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"the potential growth! that is multiplied by some growth reducer functions and a belowground competition function, these processes are included in the main function Growth.growth!\nLeaf senescence\nAgricultural defoliation","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.growth!","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.growth!","page":"Growth","title":"RegionalGrasslandSim.Growth.growth!","text":"growth!(; t, p, calc, biomass, WR)\n\nCalculates the actual growth of the plant species.\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"","category":"page"},{"location":"Modelling_API/Growth/#pot_growth","page":"Growth","title":"Potential growth","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.potential_growth!\nGrowth.calculate_LAI","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.potential_growth!","page":"Growth","title":"RegionalGrasslandSim.Growth.potential_growth!","text":"potential_growth!(; calc, SLA, biomass, PAR, potgrowth_included)\n\nCalculates the potential growth of all plant species in a specific patch.\n\nThis function is called each time step (day) for each patch. The PAR value is the photosynthetically active radiation of the day.\n\nFirst, the leaf area indices of all species are calculated (see Growth.calculate_LAI). Then, the total leaf area is computed. An inverse exponential function is used to calculate the total primary production:\n\ntexttotalgrowth =\n    PAR cdot RUE_max cdot (1 -  textexp(-alpha cdot textLAItot))\n\nThis primary production is then multiplied with the share of the leaf area index of the individual species\n\n(Image: Influence of the specific leaf area on the potential growth)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.calculate_LAI","page":"Growth","title":"RegionalGrasslandSim.Growth.calculate_LAI","text":"calculate_LAI(; SLA, biomass, LAIs)\n\nCalculate the leaf area index of all species of one habitat patch.\n\nbeginalign\ntextLAI = textSLA cdot textbiomass cdot textLAM \ntextLAI_texttot = sum textLAI\nendalign\n\nSLA specific leaf area [m² g⁻¹]\nLAM Proportion of laminae in green biomass [unitless], the value 0.62 is derived by [1]\nbiomass [kg ha⁻¹]\n\nThere is a unit conversion from the SLA and the biomass to the unitless LAI involved.\n\nThe array LAIs is mutated inplace.\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"","category":"page"},{"location":"Modelling_API/Growth/#reducer_functions","page":"Growth","title":"Reducer functions","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"The growth of each plant species in each patch is dependent on... ","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"☀ the photosynthetically active radiation Growth.radiation_reduction\nthe height of the plants in relation to the community weighted mean height Growth.height_influence!\n🌡 the air temperature Growth.temperature_reduction\n💧 the soil water content\nthe plant-available nutrients\n📈 a seasonal effect, that is modelled by the accumulated degree days Growth.seasonal_reduction","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.radiation_reduction\nGrowth.height_influence!\nGrowth.community_weighted_mean_height\nGrowth.temperature_reduction\nGrowth.water_reduction!\nGrowth.nutrient_reduction!\nGrowth.seasonal_reduction","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.radiation_reduction","page":"Growth","title":"RegionalGrasslandSim.Growth.radiation_reduction","text":"radiation_reduction(; PAR, radiation_red)\n\nReduction of radiation use efficiency at light intensities higher than 5 MJcdot m^-2cdot d^-1\n\ntextRred = textmin(1 1-gamma_1(textPAR(t) - gamma_2))\n\nThe equations and the parameter values are taken from [2].\n\nγ₁ is the empirical parameter for a decrease in RUE for high PAR values, here set to 0.0445 [m² d MJ⁻¹]\nγ₂ is the threshold value of PAR from which starts a linear decrease in RUE, here set to 5 [MJ m⁻² d⁻¹]\n\ncomment to the equation/figure: PAR values are usually between 0 and 15 MJcdot m^-2cdot d^-1 and therefore negative values of Rred are very unlikely (Image: Image of the radiation reducer function)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.height_influence!","page":"Growth","title":"RegionalGrasslandSim.Growth.height_influence!","text":"height_influence!(;\n    calc, biomass, height, height_included, height_strength = 0.5)\n\ntextheightinfluence =\n    1 +\n    fractextheightcdottextheight_textstrengthtextheight_textcwm\n    -textheight_textstrength\n\nheight_strength lies between 0 (no influence) and 1 (strong influence of the plant height)\nthe community weighted mean height height_cwm is calculated by Growth.community_weighted_mean_height\n\nIn these plots all three plant species have an equal biomass: (Image: ) (Image: )\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.community_weighted_mean_height","page":"Growth","title":"RegionalGrasslandSim.Growth.community_weighted_mean_height","text":"community_weighted_mean_height(; calc, biomass, height)\n\ntextheight_textcwm =\n    fracsum textbiomass\n    cdot textheightsum textbiomass\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.temperature_reduction","page":"Growth","title":"RegionalGrasslandSim.Growth.temperature_reduction","text":"temperature_reduction(; T, temperature_red)\n\nReduction of the potential growth if the temperature is low or too high with a step function.\n\ntexttemperature_reduction(T) =\n    begincases\n    0  textif  T  T_0 \n    fracT - T_0T_1 - T_0  textif  T_0  T  T_1 \n    1  textif  T_1  T  T_2 \n    fracT_3 - TT_3 - T_2  textif  T_2  T  T_3 \n    0  textif  T  T_3 \n    endcases\n\nEquations are taken from [3] and theses are based on [2]. T₁ is in [3] a species specific parameter, but here it is set to 12°C for all species.\n\nT₀ is the lower temperature threshold for growth, here set to 3°C\nT₁ is the lower bound for the optimal temperature for growth, here set to 12°C\nT₂ is the upper bound for the optiomal temperature for growth, here set to 20°C\nT₃ is the maximum temperature for growth, here set to 35°C\n\n(Image: Image of the temperature reducer function)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.water_reduction!","page":"Growth","title":"RegionalGrasslandSim.Growth.water_reduction!","text":"water_reduction!(;\n    calc,\n    fun_response,\n    WR,\n    water_red,\n    PET,\n    PWP,\n    WHC)\n\nSee for details: Water stress\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.nutrient_reduction!","page":"Growth","title":"RegionalGrasslandSim.Growth.nutrient_reduction!","text":"nutrient_reduction!(;\n    calc,\n    fun_response,\n    nutrient_red,\n    nutrients)\n\nSee for details: Nutrient stress\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.seasonal_reduction","page":"Growth","title":"RegionalGrasslandSim.Growth.seasonal_reduction","text":"seasonal_reduction(; ST, season_red)\n\nReduction of growth due to seasonal effects. The function is based on the yearly cumulative sum of the daily mean temperatures (ST).\n\ntextseasonal(ST) =\n    begincases\n    SEA_min  textif  ST  200 \n    SEAₘᵢₙ + (SEAₘₐₓ - SEAₘᵢₙ) * fracST - 200ST₁ - 400 \n        textif  200  ST  ST₁ - 200 \n    SEA_max  textif  ST₁ - 200  ST  ST₁ - 100 \n    SEAₘᵢₙ + (SEAₘᵢₙ - SEAₘₐₓ) * fracST - ST₂ST₂ - ST₁ - 100 \n        textif  ST₁ - 100  ST  ST₂ \n    SEA_min  textif  ST  ST₂ \n    endcases\n\nThis empirical function was developed by [1]. In contrast to [1] SEAₘᵢₙ, SEAₘₐₓ, ST₁ and ST₂ are not species specific parameters, but are fixed for all species. The values of the parameters are based on [1] and were chosen to resemble the mean of all functional groups that were described there.\n\nA seasonal factor greater than one means that growth is increased by the use of already stored resources. A seasonal factor below one means that growth is reduced as the plant stores resources [1].\n\nST is the yearly cumulative sum of the daily mean temperatures\nSEAₘᵢₙ is the minimum value of the seasonal effect, here set to 0.67 [-]\nSEAₘₐₓ is the maximum value of the seasonal effect, here set to 1.33 [-]\nST₁ and ST₂ are parameters that describe the thresholds of the step function,  here set to 625 and 1300 [°C d]\n\n(Image: Image of the seasonal effect function)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"–","category":"page"},{"location":"Modelling_API/Growth/#below_competition","page":"Growth","title":"Below-ground competition","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.below_ground_competition!","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.below_ground_competition!","page":"Growth","title":"RegionalGrasslandSim.Growth.below_ground_competition!","text":"below_ground_competition!(;\n    below,\n    traitsimilarity_biomass,\n    biomass, below_included,\n    trait_similarity,\n    below_competition_strength)\n\nModels the below-ground competiton between plant.\n\nPlant growth is reduced if a large biomass of plant species with similar functional traits is already present. The below_competition factor has a value between 0 and 1. For plant species i with N plant species present it is defined as follows:\n\ntextbelow_competition_i =\n    expleft(-fractextbelow_competition_strength1000 cdot\n        leftsum_u=1^u=N texttrait_similarity_iu cdot textbiomass_uright\n    right)\n\nThe below_competition_strength can therefore be seen as a parameter that controls the density dependence.\n\nThe trait_similarity is computed before the start of the simulation (calculation of traits similarity). and includes the traits arbuscular mycorrhizal colonisation rate (AMC) and the root surface area devided by the above ground biomass (SRSA_above).\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"","category":"page"},{"location":"Modelling_API/Growth/#Leaf-senescence","page":"Growth","title":"Leaf senescence","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.senescence!\nGrowth.seasonal_component_senescence","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.senescence!","page":"Growth","title":"RegionalGrasslandSim.Growth.senescence!","text":"senescence!(; sen, ST, biomass, μ)\n\nbeginalign\n    LL = 10 ^  left(log10(SLA)  - 241right)  -038 cdotfrac3652512 \n    μ = fractextsen_intercept1000 + fractextsen_rate1000 cdot frac1LL \n    textsenescence = μ cdot textSEN cdot textbiomass\nendalign\n\nLL leaf life span [d]\nSLA specific leaf area [fraccm^2g] rightarrow this includes a unit conversion of the SLA values (in the model the unit of SLA is fracm^2g)\nμ leaf senescence rate [frac1d]\nSEN seasonal component of the senescence (between 1 and 3)\nsen_intercept α value of a linear equation that models the influence of the leaf senescence rate μ on the total senescence rate\nsen_rate β value of a linear equation that models the influence of the leaf senescence rate μ on the total senescence rate\n\nThe parameters textsen_intercept and textsen_rate were divided by 1000 to avoid very low numbers.\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.seasonal_component_senescence","page":"Growth","title":"RegionalGrasslandSim.Growth.seasonal_component_senescence","text":"seasonal_component_senescence(;\n    ST,\n    Ψ₁ = 775,\n    Ψ₂ = 3000,\n    SENₘᵢₙ = 1,\n    SENₘₐₓ = 3)\n\nSeasonal factor for the senescence rate.\n\nbeginalign*\nSEN =\nbegincases\nSEN_min   textif  ST  Ψ_1 \nSEN_min+(SEN_max - SEN_min) fracST - Ψ_1Ψ_2 - Ψ_1  textif Ψ_1  ST  Ψ_2 \nSEN_max   textif ST  Ψ_2\nendcases  \nendalign*\n\nST yearly accumulated degree days [C]\nΨ₁=775  [Ccdot d]\nΨ₂=3000 [Ccdot d]\nSEN_min=1\nSEN_max=3\n\n(Image: Seasonal component death rate)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"","category":"page"},{"location":"Modelling_API/Growth/#Agricultural-defoliation","page":"Growth","title":"Agricultural defoliation","text":"","category":"section"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Biomass is removed by...","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"🐄 Growth.grazing! and Growth.trampling!\n🚜 Growth.mowing!","category":"page"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"Growth.grazing!\nGrowth.mowing!\nGrowth.trampling!","category":"page"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.grazing!","page":"Growth","title":"RegionalGrasslandSim.Growth.grazing!","text":"grazing!(; calc, LD, biomass, ρ, grazing_half_factor)\n\nbeginalign\nμₘₐₓ = κ cdot textLD \nh = frac1μₘₐₓ \na = frac1textgrazing_half_factor^2 cdot h \ntexttotgraz = fraca cdot (sum textbiomass)^2\n                    1 + acdot hcdot (sum textbiomass)^2 \ntextgraz = texttotgraz cdot fracρ cdot textbiomass\n                    sum ρ cdot textbiomass\nendalign\n\nLD daily livestock density (livestock units ha⁻¹)\nκ daily consumption of one livestock unit, follows [4]\nρ appetence of the plant species for livestock, dependent on nitrogen per leaf mass (LNCM)\ngrazing_half_factor is the half-saturation constant\nequation partly based on [3]\n\nInfluence of grazing (livestock density = 2), all plant species have an equal amount of biomass (total biomass / 3): (Image: Image of grazing effect)\n\nInfluence of grazing_half_factor (LD is set to 2): (Image: )\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.mowing!","page":"Growth","title":"RegionalGrasslandSim.Growth.mowing!","text":"mowing!(;\n    calc,\n    mowing_height,\n    days_since_last_mowing,\n    height,\n    biomass,\n    mowing_mid_days)\n\nbeginalign\n    lambda = fractextmown_heightheight\n    textmow_factor = frac11+exp(-01*(textdays_since_last_mowing\n        - textmowing_mid_days)\n    textmow = lambda cdot textbiomass\nendalign\n\nThe mow_factor has been included to account for the fact that less biomass is mown when the last mowing event was not long ago. Influence of mowing for plant species with different heights (height): (Image: Image of mowing effect)\n\nVisualisation of the mow_factor: (Image: )\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/#RegionalGrasslandSim.Growth.trampling!","page":"Growth","title":"RegionalGrasslandSim.Growth.trampling!","text":"trampling!(; calc, LD, biomass, height, trampling_factor)\n\nbeginalign\nω = fractexttrampling_factorheight^025 \ntexttrampled_biomass =\n    begincases\n        textbiomass\n             textif LD  ω\n        fractextbiomass2\n        cdot left(1- cosleft(fracπcdottextLDωright)right)\n             textotherwise\n\tendcases\nendalign\n\nIt is assumed that tall plants (trait: height) are stronger affected by trampling. A cosine function is used to model the influence of trampling.\n\nIf the livestock density is higher than ω, all the biomass of that plant species will be removed. This is unlikely to be the case.\n\nbiomass [frackgha]\nLD daily livestock density [fractextlivestock unitsha]\ntrampling_factor [ha]\nheight canopy height [m]\n\n(Image: Image of trampling effect)\n\n\n\n\n\n","category":"function"},{"location":"Modelling_API/Growth/","page":"Growth","title":"Growth","text":"","category":"page"},{"location":"Modelling_API/Difference_equation/#Difference-equation","page":"Difference equation","title":"Difference equation","text":"","category":"section"},{"location":"Modelling_API/Difference_equation/","page":"Difference equation","title":"Difference equation","text":"Modules = [RegionalGrasslandSim]\nPublic = true\nPrivate = true","category":"page"},{"location":"Modelling_API/Difference_equation/#RegionalGrasslandSim.initialize_parameters-Tuple{}","page":"Difference equation","title":"RegionalGrasslandSim.initialize_parameters","text":"initialize_parameters(; input_obj)\n\nInitialize all parameters and store them in one object.\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Difference_equation/#RegionalGrasslandSim.main_loop!-Tuple{}","page":"Difference equation","title":"RegionalGrasslandSim.main_loop!","text":"main_loop!(; calc, ts, p)\n\nRun the main loop for all days.\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Difference_equation/#RegionalGrasslandSim.one_day!-Tuple{}","page":"Difference equation","title":"RegionalGrasslandSim.one_day!","text":"one_day!(; calc, p, t)\n\nCalculate the density differences of all state variables of one day.\n\n\n\n\n\n","category":"method"},{"location":"Modelling_API/Difference_equation/#RegionalGrasslandSim.solve_prob-Tuple{}","page":"Difference equation","title":"RegionalGrasslandSim.solve_prob","text":"solve_prob(; input_obj, inf_p, calc = nothing)\n\nSolve the model for one site.\n\nIntialize the parameters, the state variables and the output vectors. In addition some vectors are preallocated to avoid allocations in the main loop. Then, run the main loop and store the results in the output vectors.\n\n\n\n\n\n","category":"method"},{"location":"#RegionalGrasslandSim.jl","page":"Home","title":"RegionalGrasslandSim.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RegionalGrasslandSim.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"Depth = 5","category":"page"}]
}
