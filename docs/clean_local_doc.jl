import Pkg
Pkg.rm("RegionalGrasslandSim")
Pkg.rm("RegionalGrasslandVis")

img_path = "docs/src/img/"
[img[1] == '.' ? "" : rm("$(img_path)$img") for img in readdir(img_path)]
rm("docs/build/", recursive = true)
