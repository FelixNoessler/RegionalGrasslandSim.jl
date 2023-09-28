# include("../lib/RegionalGrasslandValid/scripts/optimization.jl")

using Metaheuristics
using RegionalGrasslandValid
import StatsBase
import RegionalGrasslandSim as sim

function run_optimization(;
    f_calls_limit = Inf,
    time_limit = Inf,
    iterations = 10,
    explos = ["HEG", "AEG", "SEG"])
    param_names = [
        "moistureconv_alpha", "moistureconv_beta",
        "senescence_intercept", "senescence_rate",
        "below_competition_strength", "trampling_factor", "grazing_half_factor",
        "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
        "max_AMC_nut_reduction", "max_SRSA_nut_reduction",
        "b_biomass",
        "b_SLA", "b_LNCM", "b_AMC", "b_height", "b_SRSA_above",
        "b_soilmoisture"]

    training_plots = ["$(explo)$(lpad(i, 2, "0"))" for i in 1:9 for explo in explos]

    function ll_batch(X)
        # selected_plots = StatsBase.sample(training_plots, batch_size; replace = false)
        selected_plots = training_plots
        N = size(X, 1)
        ll_mat = Array{Float64}(undef, N, length(selected_plots))

        Threads.@threads for i in 1:N
            inf_p = (; zip(Symbol.(param_names), X[i, :])...)

            for (p, plotID) in enumerate(selected_plots)
                ll_mat[i, p] = loglikelihood_model(sim;
                    nspecies = 25,
                    plotID,
                    inf_p,
                    pretty_print = false)
            end
        end

        ### free RAM space
        GC.gc()

        return vec(sum(abs.(ll_mat); dims = 2))
    end

    #         mc_α mc_β s_i s_r below tram graz mow SRSA SLA  AMC  SRSA_n
    lb_prep = [0, 0, 0, 0, 0, 100, 0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ub_prep = [80, 300, 10, 10, 1, 300, 2000, 50, 1.0, 1.0, 1.0, 1.0]

    nscale_params = 7
    lb_b = zeros(nscale_params)
    ub_b = fill(5e3, nscale_params)
    lb = vcat(lb_prep, lb_b)
    ub = vcat(ub_prep, ub_b)

    bounds = boxconstraints(; lb, ub)

    options = Options(;
        f_calls_limit,
        iterations,
        time_limit,
        store_convergence = true,
        debug = true,
        parallel_evaluation = true)

    live_plot(st) = begin
        @show vals = st.best_sol.x
    end

    ## DE or ECA
    algorithm = DE(; options = options)
    result = optimize(ll_batch, bounds, algorithm, logger = live_plot)

    @show minimizer(result)

    return result
end

optim_result = run_optimization(;
    iterations = 10)

# let
#     f_calls, best_f_value = convergence(optim_result)
#     fig, _ = lines(f_calls, best_f_value, label = "ECA")

#     # f_calls, best_f_value = convergence(optim_result1)
#     # lines!(f_calls, best_f_value, label="DE")

#     # f_calls, best_f_value = convergence(optim_result2)
#     # lines!(f_calls, best_f_value, label="DE_40")

#     axislegend()
#     fig
# end

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
