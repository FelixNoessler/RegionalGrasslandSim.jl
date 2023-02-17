using RegionalGrasslandSim
using Test

@testset "RegionalGrasslandSim.jl" begin
    @test RegionalGrasslandSim.greet_your_package_name() == "Hello RegionalGrasslandSim!"
    @test RegionalGrasslandSim.greet_your_package_name() != "Hello world!"
end

