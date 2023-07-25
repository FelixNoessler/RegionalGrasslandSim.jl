####### build the documentation locally
# cd RegionalGrasslandSim
# julia --project=docs/ --startup-file=no
# using Revise; import Pkg; Pkg.instantiate(); Pkg.develop(path="lib/RegionalGrasslandVis/"); Pkg.develop(path="."); include("docs/make.jl")
## to redo the documentation:
# include("docs/make.jl")
## to clean everything for commits/push:
# include("docs/clean_local_doc.jl")

using Documenter
import RegionalGrasslandSim
import RegionalGrasslandSim as sim
import RegionalGrasslandVis as vis

####### Create Bilbiography
using DocumenterCitations
bib = CitationBibliography("docs/src/lit.bib")

####### create images for the document
#### functional response
vis.potential_growth(sim; path = "docs/src/img/sla_potential_growth.svg")
vis.srsa_water_response(sim;
    path = "docs/src/img/srsa_water_response.svg",
    max_SRSA_water_reduction = 1)
vis.srsa_water_response(sim;
    path = "docs/src/img/srsa_water_response_0_5.svg",
    max_SRSA_water_reduction = 0.5)
vis.srsa_nut_response(sim;
    path = "docs/src/img/srsa_nut_response.svg",
    max_SRSA_nut_reduction = 1)
vis.srsa_nut_response(sim;
    path = "docs/src/img/srsa_nut_response_0_5.svg",
    max_SRSA_nut_reduction = 0.5)
vis.amc_nut_response(sim;
    path = "docs/src/img/amc_nut_response.svg",
    max_AMC_nut_reduction = 1)
vis.amc_nut_response(sim;
    path = "docs/src/img/amc_nut_response_0_5.svg",
    max_AMC_nut_reduction = 0.5)
vis.sla_water_response(sim;
    path = "docs/src/img/sla_water_response.svg",
    max_SLA_water_reduction = 1.0)
vis.sla_water_response(sim;
    path = "docs/src/img/sla_water_response_0_5.svg",
    max_SLA_water_reduction = 0.5)

#### reducer functions
vis.temperatur_reducer(sim; path = "docs/src/img/temperature_reducer.svg")
vis.radiation_reducer(sim; path = "docs/src/img/radiation_reducer.svg")

#### seasonal effects
vis.seasonal_effect(sim; path = "docs/src/img/seasonal_reducer.svg")
vis.seasonal_component_senescence(sim;
    path = "docs/src/img/seasonal_factor_senescence.svg")

#### land use
vis.mowing(sim; path = "docs/src/img/mowing.svg")
vis.mow_factor(; path = "docs/src/img/mow_factor.svg")
vis.grazing(sim; path = "docs/src/img/grazing.svg")
vis.grazing_half_factor(; path = "docs/src/img/grazing_half_factor.svg")
vis.trampling(sim; path = "docs/src/img/trampling.svg")

# for prettyurls you need locally a live server
makedocs(bib;
    sitename = "RegionalGrasslandSim.jl",
    format = Documenter.HTML(prettyurls = true,
        mathengine = MathJax3()),
    modules = [RegionalGrasslandSim],
    pages = Any["Home" => "index.md",
        "Modelling API" => Any["Difference equation" => "Modelling_API/Difference_equation/index.md",
            "Growth" => "Modelling_API/Growth/index.md",
            "Water dynamics" => "Modelling_API/Water_dynamics/index.md",
            "Functional response" => "Modelling_API/Functional_response/index.md"],
        "References" => "References.md"])

deploydocs(repo = "github.com/FelixNoessler/RegionalGrasslandSim.jl.git")
