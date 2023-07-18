import Pkg
Pkg.rm("RegionalGrasslandSim")
Pkg.rm("RegionalGrasslandVis")
run(`rm docs/src/img/\*.svg`)
run(`rm -r docs/build/`)
