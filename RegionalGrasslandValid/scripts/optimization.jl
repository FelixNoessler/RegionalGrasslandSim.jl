using Metaheuristics
import RegionalGrasslandSim as sim
import RegionalGrasslandValid as valid



function run_optimization(; f_calls_limit=1000, time_limit=200.0)
    function f(x)
        param_names = [
            "sigma_biomass","sigma_evaporation","sigma_soilmoisture",
            "moisture_conv","senescence_intercept","senescence_rate",
            "below_competition_strength","trampling_factor","grazing_half_factor",
            "mowing_mid_days","max_SRSA_water_reduction","max_SLA_water_reduction",
            "max_AMC_nut_reduction","max_SRSA_nut_reduction"]

        inf_p = (; zip(Symbol.(param_names), x)...)
        ll = loglikelihood_model(sim;
            plotID="AEG02", inf_p)
        # ll = ll_VIPS(sim; inf_p)

        return abs(ll)
    end

    l_bounds = [0,    0,   0,   0.1, 1e-6, 0,   0.0,   50,  250, 10, 0.0, 0.0, 0.0, 0.0]
    u_bounds = [1e10, 100, 100, 1.2, 1,    100, 0.1,   200, 550, 25, 1.0, 1.0, 1.0, 1.0]
    bounds = [l_bounds u_bounds]'
    options = Options(; f_calls_limit, time_limit)

    ## DE or ECA
    algorithm = ECA(options = options)
    result = optimize(f, bounds, algorithm)

    @show minimizer(result)

    return result
end

run_optimization()
