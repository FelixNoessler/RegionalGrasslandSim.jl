module RegionalGrasslandSim

using JLD2
using Unitful
using Distributions
import Dates
import NeutralLandscapes

export solve_prob

include("neighbours.jl")
include("input_prep.jl")
include("Growth/Growth.jl")
include("Water/Water.jl")
include("Functional response/FunctionalResponse.jl")
include("main_functions.jl")
include("Traits/Traits.jl")

function __init__()
    datapath = joinpath(@__DIR__, "..", "lib", "RegionalGrasslandData")
    Traits.load_gm(datapath)
end

end
