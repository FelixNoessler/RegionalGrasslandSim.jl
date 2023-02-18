using Documenter

push!(LOAD_PATH,"../src/")
using RegionalGrasslandSim

makedocs(
    sitename = "RegionalGrasslandSim",
    format = Documenter.HTML(),
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
