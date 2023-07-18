library(BayesianTools)

out <- readRDS("../RegionalGrasslandValid/saved_chains/2023_07_15_dream_14_1000_HEG01.rds")

# out
# getSample(out)
summary(out)


tracePlot(out)
# correlationPlot(out)
marginalPlot(out, prior = TRUE, start=1)
# marginalLikelihood(out)



cat(MAP(out)$parametersMAP, sep=", ")
