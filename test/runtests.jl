using RegionalGrasslandSim
using Test

@testset "RegionalGrasslandSim.jl" begin
    @test RegionalGrasslandSim.Growth.greet_your_package_name() == "Hello RegionalGrasslandSim!"
    @test RegionalGrasslandSim.Growth.greet_your_package_name() != "Hello world!"
end

