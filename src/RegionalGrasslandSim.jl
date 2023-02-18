module RegionalGrasslandSim

using OrdinaryDiffEq
import ComponentArrays as ca

export greet_your_package_name

export discrete_prob
include("Growth/Growth.jl")
include("main_functions.jl")

end
