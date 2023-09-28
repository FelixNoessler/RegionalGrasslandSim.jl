library(JuliaCall)
library(BayesianTools)
library(readr)

prior_df <- read_csv("lib/RegionalGrasslandData/validation/priors.csv")


prepare_julia <- function(){
    julia_command('using Revise')
    julia_command('if isdir("Grassland_data_analysis") cd("Grassland_data_analysis"); end')
    julia_command('import Pkg; Pkg.activate(".")')
    julia_command('import RegionalGrasslandSim as sim')
    julia_command('using RegionalGrasslandValid')
}

#--------------------------------------------
parameter_names <- prior_df$param_name
prior <- createTruncatedNormalPrior(
    mean = prior_df$mean,
    sd = prior_df$sd,
    lower = prior_df$lb, 
    upper = prior_df$ub)

ll_VIPs <- function(param){
    param_list <- split(param, parameter_names)
    loglik <- julia_call(
            "ll_VIPS_t", 
            julia_eval("sim"),
            inf_p = param_list,
            nspecies = 25)
    return(loglik)
}

ll <- function(param){
    param_list <- split(param, parameter_names)
    loglik <- julia_call(
            "loglikelihood_model", 
            julia_eval("sim"),
            plotID = "HEG01",
            inf_p = param_list,
            nspecies = 25)
    return(loglik)
}


prepare_julia()

# julia_command("GC.gc()")
# test_p <- prior$sampler()
# start_time <- Sys.time()
# ll(test_p)
# end_time <- Sys.time()
# end_time - start_time


bayesianSetup = createBayesianSetup(
    likelihood = ll_VIPs, 
    prior = prior,
    names = parameter_names,
    parallel = F)

#############
out <- runMCMC(
    bayesianSetup = bayesianSetup, 
    sampler = "DEzs",  # DREAMzs , default: DEzs
    settings = list(
        iterations = 10000, 
        message = T, 
        nrChains = 3,
        burnin = 0))
plot(out)
filename = paste0(
    "../tmp/" , 
    format(Sys.time(), "%Y_%m_%d"), 
    length(parameter_names),
    ".rds")
saveRDS(out, file = filename)

