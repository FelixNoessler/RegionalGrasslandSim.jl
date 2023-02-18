using Documenter
using RegionalGrasslandSim

# for prettyurls you need locally a live server
makedocs(
    sitename = "RegionalGrasslandSim",
    format = Documenter.HTML(prettyurls = true), 
    modules = [RegionalGrasslandSim],
    pages = Any[
        "Home" => "index.md",
        "Modelling API"=> Any[
            "Diffrence equation" => "Modelling API/Difference equation/index.md",
            "Growth" => "Modelling API/Growth/index.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/FelixNoessler/RegionalGrasslandSim.jl.git",
)
