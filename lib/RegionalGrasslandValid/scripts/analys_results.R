library(BayesianTools)

getwd()
out <- readRDS("tmp/2023_09_2219.rds")
out_part <- getSample(out, start = 2000)
# out
# getSample(out)
summary(out)


traceplot(out)
correlationPlot(out)
marginalPlot(out, prior = TRUE, start=2000)


getSample(out, 1)

cat(MAP(out)$parametersMAP, sep=", ")
