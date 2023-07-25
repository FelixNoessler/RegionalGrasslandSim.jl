library(JuliaCall)
library(BayesianTools)

prepare_julia <- function(){
    # julia_command('if isdir("Grassland_data_analysis") cd("Grassland_data_analysis"); end')
    # julia_command('import Pkg; Pkg.activate(".")')
    julia_command('import RegionalGrasslandSim as sim')
    julia_command('using RegionalGrasslandValid')
}

#--------------------------------------------
parameter_names <- c(
    "sigma_biomass",
    "sigma_evaporation",
    "sigma_soilmoisture",
    "moisture_conv",
    "senescence_intercept",
    "senescence_rate",
    "below_competition_strength",
    "trampling_factor", 
    "grazing_half_factor",
    "mowing_mid_days",
    "max_SRSA_water_reduction",
    "max_SLA_water_reduction",
    "max_AMC_nut_reduction",
    "max_SRSA_nut_reduction")

prior <- createTruncatedNormalPrior(
            # σ_bio σ_ev σ_mo m_c  s_i   s_r  below  tram graz  mow  SRSA SLA  AMC  SRSA_n
    mean =  c(0,    0,   0,  0.8, 0.001, 5,   0.001, 100, 1500, 40,  0.5, 0.5, 0.5, 0.5),
    sd =    c(1e4,  10,  10,  0.2, 0.1,  2,   0.01,  50,  500,  40,  0.2, 0.2, 0.2, 0.2),
    lower = c(0,    0,   0,   0.1, 1e-6, 0,   0.0,   50,  500,  10,  0.0, 0.0, 0.0, 0.0), 
    upper = c(1e10, 100, 100, 1.2, 1,    100, 0.1,   200, 5000, 150, 1.0, 1.0, 1.0, 1.0))

ll_VIPs <- function(param){
    param_list <- split(param, parameter_names)
    loglik <- julia_call(
            "ll_VIPS_t", 
            julia_eval("sim"),
            inf_p=param_list)
    return(loglik)
}

ll <- function(param){
    param_list <- split(param, parameter_names)
    loglik <- julia_call(
            "loglikelihood_model", 
            julia_eval("sim"),
            plotID="HEG01",
            inf_p=param_list)
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
    likelihood = ll, 
    prior = prior,
    names = parameter_names,
    parallel = F)

############# DREAMzs , DEzs
out <- runMCMC(
    bayesianSetup = bayesianSetup, 
    sampler = "DREAMzs",  # DREAMzs , DEzs
    settings = list(
        iterations = 3000, 
        message = T, 
        nrChains = 1,
        burnin = 0))

filename = paste0(
    "/home/felix/Dokumente/Promotion/RegionalGrasslandValid/saved_chains/" , 
    format(Sys.time(), "%Y_%m_%d_dream_"), 
    length(parameter_names),
    "_1000_HEG01.rds")
saveRDS(out, file = filename)

