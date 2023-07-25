using Metaheuristics
using RegionalGrasslandValid
import StatsBase
import RegionalGrasslandSim as sim

function run_optimization(;
    f_calls_limit = Inf,
    time_limit = Inf,
    iterations = 10,
    selected_plots = NaN,
    batch_size,
    explos = ["HAI", "AEG", "SCH"])
    param_names = [
        "sigma_biomass", "sigma_evaporation", "sigma_soilmoisture",
        "moisture_conv", "senescence_intercept", "senescence_rate",
        "below_competition_strength", "trampling_factor", "grazing_half_factor",
        "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
        "max_AMC_nut_reduction", "max_SRSA_nut_reduction",
    ]

    training_plots = ["$(explo)$(lpad(i, 2, "0"))" for i in 1:50 for explo in explos]
    function ll_batch(X)
        selected_plots = StatsBase.sample(training_plots, batch_size; replace = false)
        @info selected_plots
        N = size(X, 1)
        ll_mat = Array{Float64}(undef, N, length(selected_plots))

        Threads.@threads for i in 1:N
            inf_p = (; zip(Symbol.(param_names), X[i, :])...)
            for (p, plotID) in enumerate(selected_plots)
                ll_mat[i, p] = loglikelihood_model(sim;
                    plotID,
                    inf_p)
            end
        end

        ### free RAM space
        GC.gc() #ccall(:malloc_trim, Cvoid, (Cint,), 0)

        return vec(sum(abs.(ll_mat); dims = 2))
    end

    #    σ_bio  σ_ev σ_mo m_c  s_i s_r  below tram graz mow SRSA SLA  AMC  SRSA_n
    lb = [1000, 0, 0, 0.1, 0, 0, 0.0, 50, 250, 5, 0.0, 0.0, 0.0, 0.0]
    ub = [1010, 100, 100, 1.2, 5, 100, 15, 200, 1000, 25, 1.0, 1.0, 1.0, 1.0]
    bounds = boxconstraints(; lb, ub)

    options = Options(;
        f_calls_limit,
        iterations,
        time_limit,
        store_convergence = true,
        debug = true,
        parallel_evaluation = true)

    ## DE or ECA
    D = length(param_names)
    K = 7
    algorithm = ECA(; η_max = 2.0, K, N = 100, adaptive = true, options = options)
    result = optimize(ll_batch, bounds, algorithm)
    @show minimizer(result)

    return result
end

optim_result = run_optimization(;
    explos = ["HEG"],
    batch_size = 2,
    iterations = 25)

using GLMakie
let
    f_calls, best_f_value = convergence(optim_result)
    fig, _ = lines(f_calls, maximum.(best_f_value), label = "ECA")

    # f_calls, best_f_value = convergence(optim_result1)
    # lines!(f_calls, best_f_value, label="DE")

    # f_calls, best_f_value = convergence(optim_result2)
    # lines!(f_calls, best_f_value, label="DE_40")

    axislegend()
    fig
end

############ playgorund
############
############

# function f_parallel_single(X)
#     ## set parallel_evaluation to true!
#     N = size(X, 1)
#     lls = Array{Float64}(undef, N)
#     Threads.@threads for i in 1:N
#         inf_p = (; zip(Symbol.(param_names), X[i, :])...)

#         ll = loglikelihood_model(sim;
#             plotID=rand(selected_plots),
#             inf_p)
#         lls[i] = ll
#     end
#     ### free RAM space
#     GC.gc()
#     # ccall(:malloc_trim, Cvoid, (Cint,), 0)

#     return abs.(lls)
# end
# training_plots = ["$(explo)$(lpad(i, 2, "0"))" for i in 1:50 for explo in explos]
# function ll_batch(X)
#     selected_plots = rand(training_plots, batch_size)

#     N = size(X, 1)
#     ll_mat = Array{Float64}(undef, N, length(selected_plots))

#     # lk = Threads.SpinLock()
#     Threads.@threads for i in 1:N
#         inf_p = (; zip(Symbol.(param_names), X[i, :])...)
#         for (p,plotID) in enumerate(selected_plots)
#             ll_mat[i,p] = loglikelihood_model(sim;
#                 plotID,
#                 inf_p)
#             # lock(lk)
#             # try
#             #     ll_mat[i,p] = ll
#             # catch e
#             #     println(e)
#             # finally
#             #     unlock(lk)
#             # end
#         end
#     end

#     ### free RAM space
#     GC.gc()
#     # ccall(:malloc_trim, Cvoid, (Cint,), 0)

#     return vec(sum(abs.(ll_mat); dims=2))
# end
# function f_parallel_single(X)
#     ## set parallel_evaluation to true!
#     N = size(X, 1)
#     lls = Array{Float64}(undef, N)
#     Threads.@threads for i in 1:N
#         inf_p = (; zip(Symbol.(param_names), X[i, :])...)

#         ll = loglikelihood_model(sim;
#             plotID=rand(selected_plots),
#             inf_p)
#         lls[i] = ll
#     end
#     ### free RAM space
#     GC.gc()
#     # ccall(:malloc_trim, Cvoid, (Cint,), 0)

#     return abs.(lls)
# end

# function f_parallel_all(X)
#     ## set parallel_evaluation to true!
#     N = size(X, 1)
#     ll_mat = Array{Float64}(undef, N, length(VIP_plots))

#     Threads.@threads for i in 1:N
#         inf_p = (; zip(Symbol.(param_names), X[i, :])...)
#         for (p,plotID) in enumerate(VIP_plots)
#             ll_mat[i,p] = loglikelihood_model(sim;
#                 plotID,
#                 inf_p)
#         end
#     end

#     ### free RAM space
#     GC.gc()
#     ccall(:malloc_trim, Cvoid, (Cint,), 0)

#     return vec(sum(abs.(ll_mat); dims=2))
# end
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
