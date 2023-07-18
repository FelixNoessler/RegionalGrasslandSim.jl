using Metaheuristics
using RegionalGrasslandValid
import RegionalGrasslandSim as sim


function run_optimization(;
    f_calls_limit=Inf,
    time_limit=Inf,
    iterations=10,
    plotID=nothing)

    param_names = [
        "sigma_biomass","sigma_evaporation","sigma_soilmoisture",
        "moisture_conv","senescence_intercept","senescence_rate",
        "below_competition_strength","trampling_factor","grazing_half_factor",
        "mowing_mid_days","max_SRSA_water_reduction","max_SLA_water_reduction",
        "max_AMC_nut_reduction","max_SRSA_nut_reduction"]

    # function f_parallel_single(X)
    #     ## set parallel_evaluation to true!
    #     N = size(X, 1)
    #     lls = Array{Float64}(undef, N)
    #     Threads.@threads for i in 1:N
    #         inf_p = (; zip(Symbol.(param_names), X[i, :])...)

    #         ll = loglikelihood_model(sim;
    #             plotID,
    #             inf_p)
    #         lls[i] = ll
    #     end
    #     ### free RAM space
    #     GC.gc()
    #     ccall(:malloc_trim, Cvoid, (Cint,), 0)

    #     return abs.(lls)
    # end

    VIP_plots = ["$(explo)0$i" for i in 1:9 for explo in ["HEG", "SEG", "AEG"]];
    function f_parallel_all(X)
        ## set parallel_evaluation to true!
        N = size(X, 1)
        ll_mat = Array{Float64}(undef, N, length(VIP_plots))

        Threads.@threads for i in 1:N
            inf_p = (; zip(Symbol.(param_names), X[i, :])...)
            for (p,plotID) in enumerate(VIP_plots)
                ll_mat[i,p] = loglikelihood_model(sim;
                    plotID,
                    inf_p)
            end
        end

        ### free RAM space
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)

        return vec(sum(abs.(ll_mat); dims=2))
    end
    # function f(x)
    #     ## set parallel_evaluation to false!
    #     inf_p = (; zip(Symbol.(param_names), x)...)

    #     ll = 0.0
    #     if isnothing(plotID)
    #         ## f_parallel_all is faster!
    #         ll = ll_VIPS_t(sim; inf_p)
    #     else
    #         ll = loglikelihood_model(sim;
    #             plotID, inf_p)
    #     end
    #     return abs(ll)
    # end

    l_bounds = [0,    0,   0,   0.1, 1e-6, 0,   0.0,   50,  250, 10, 0.0, 0.0, 0.0, 0.0]
    u_bounds = [1e10, 100, 100, 1.2, 1,    100, 0.1,   200, 550, 25, 1.0, 1.0, 1.0, 1.0]
    bounds = [l_bounds u_bounds]'

    options = Options(;
        f_calls_limit,
        iterations,
        time_limit,
        store_convergence = true,
        debug = true,
        parallel_evaluation = true)

    ## DE or ECA
    algorithm = ECA(options = options)
    result = optimize(f_parallel_all, bounds, algorithm)
    @show minimizer(result)

    return result
end

optim_result1 = run_optimization(;
    plotID=nothing,
    iterations=20,
    time_limit=60*40.0)

using CairoMakie
let
    f_calls, best_f_value = convergence(optim_result1)
    fig, _ = lines(f_calls, best_f_value, label="ECA")

    # f_calls, best_f_value = convergence(optim_result1)
    # lines!(f_calls, best_f_value, label="DE")

    # f_calls, best_f_value = convergence(optim_result2)
    # lines!(f_calls, best_f_value, label="DE_40")

    axislegend()
    fig
end
