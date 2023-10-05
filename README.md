# RegionalGrasslandSim.jl
[![][docs-dev-img]][docs-dev-url] [![SciML Code Style](https://img.shields.io/badge/code%20style-SciML-blue)](https://github.com/SciML/SciMLStyle)

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://felixnoessler.github.io/RegionalGrasslandSim.jl/dev/




```mermaid
stateDiagram-v2
    direction LR

    classDef stateStyle fill:orange ,stroke-width:2px,stroke:black;

    classDef biomassStyle fontstyle:bold, fill:DarkSeaGreen;
    classDef waterStyle fill:LightCyan;

    biomass_dyn:::biomassStyle
    water_dyn:::waterStyle

    %%%%%%%%%%%%%%%% input
    input: input creation

    %%%% climate
    climate: Climate time series
    temp: Meain air temperature [°C]
    prec: Precipitation [mm]
    par: Photosynthetically active radiation [W/m2]

    %%%% site
    site: Site properties
    nut: Nutrient availability 
    whc: Water holding capacity [mm]
    pwp: Permanent wilting point [mm]

    %%%% species
    species: Species traits
    sla: Specific leaf area [m2/kg]
    srsa: Specific root surface area per aboveground biomass [m2/kg]
    amc: Arbuscular mycorrhizal colonization [%]
    height: Plant height [m]

    %%%%%%%%%%%%%%%% main loop
    daily: loop every day

    %%%% Biomass dynamics
    biomass_dyn: Biomass dynamics
    biomass_state: State biomass [kg/ha]

    %% Growth
    pot_growth: Potential growth 
    lai: Individual leaf area index
    lai_sum: Total leaf area index

    actual_growth: Actual growth
    growth_red: Growth reduction factor
    water_red: Water reduction factor
    nut_red: Nutrient reduction factor
    radiation_red: Radiation reduction factor
    temp_red: Temperature reduction factor

    %%%% Water dynamics
    water_dyn: water dynamics
    water_state: State water [mm]
    av_water: Plant available water index


    water_state:::stateStyle
    biomass_state:::stateStyle

    state input {
        direction LR

        state climate {
            direction LR
            temp
            prec
            par
        }

        state management {
            direction LR
            grazing
            mowing
        }

        state site {
            direction LR
            nut whc pwp
        }

        state species {
            direction LR
            sla srsa amc height
        }
    }

    state daily {
        direction LR
        lai
        lai_sum
        av_water

        state biomass_dyn {
            direction LR

            biomass_state
            pot_growth
            growth_red
            water_red
            radiation_red
            temp_red
            nut_red
            height_influence
            actual_growth
        }
        
        state water_dyn {
            direction LR
            water_state
            evaporation
            transpiration
            drain
        }

    }

    %%%% water related
    lai_sum --> transpiration
    whc --> drain
    whc --> evaporation
    whc --> transpiration
    pwp --> transpiration
    prec --> water_state
    evaporation --> water_state
    transpiration --> water_state
    drain --> water_state
    water_state --> av_water

    %%%% leaf area index
    biomass_state --> lai
    lai --> lai_sum
    sla --> lai

    %%%% Growth
    actual_growth --> biomass_state
    growth_red --> actual_growth

    %% potential growth 
    lai_sum --> pot_growth
    lai --> pot_growth
    poten_growth --> actual_growth
    par --> pot_growth

    %% growth reduction due to water
    av_water --> water_red
    water_red --> growth_red
    sla --> water_red
    srsa --> water_red
    amc --> water_red

    %% growth reduction due to nutrients
    nut_red --> growth_red
    nut --> nut_red
    srsa --> nut_red
    amc --> nut_red

    %% growth factor due to height
    height --> height_influence
    biomass_state --> height_influence
    height_influence --> growth_red

    %% growth reduction due to temperature
    temp --> temp_red
    temp_red --> growth_red

    %% growth reduction due to radiation
    par --> radiation_red
    radiation_red --> growth_red









```

![model overview](assets/screenshot.png)

[ECEM 2023 slides](assets/ECEM_2023_presentation.pdf)

