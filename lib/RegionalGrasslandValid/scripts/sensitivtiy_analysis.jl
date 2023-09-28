import RegionalGrasslandSim as sim
using RegionalGrasslandValid
using RegionalGrasslandVis ## for plotting theme
using Statistics
using GLMakie
Makie.inline!(true)

########################### input preparataion
param_vals = [
    1155.299246041266,
    2.1787199839476785, 11.452466470742415,
    0.3801072247902645, 150.1232967984496, 958.214969168516,
    17.821776550046422, 0.9416295869463482, 0.9290589640747696,
    0.6261919010370469, 0.708332191866904]
param_names = [
    "sigma_biomass",
    "senescence_intercept", "senescence_rate",
    "below_competition_strength", "trampling_factor", "grazing_half_factor",
    "mowing_mid_days", "max_SRSA_water_reduction", "max_SLA_water_reduction",
    "max_AMC_nut_reduction", "max_SRSA_nut_reduction"];
inf_p = (; zip(Symbol.(param_names), param_vals)...)

function ll_batch(X; seed)
    explos = ["HEG", "AEG", "SEG"]
    selected_plots = ["$(explo)$(lpad(i, 2, "0"))" for explo in explos for i in 1:9]

    N = size(X, 2)
    ll_mat = Array{Float64}(undef, N, length(selected_plots))

    Threads.@threads for i in 1:N
        inf_p = (; zip(Symbol.(param_names), X[:, i])...)
        for (p, plotID) in enumerate(selected_plots)
            s = rand(1:10000)
            ll_mat[i, p] = loglikelihood_model(sim;
                nspecies = 25,
                plotID,
                inf_p,
                only_likelihood = true,
                seed = s)
        end
    end

    return ll_mat
end

let
    seed = 100
    nvals = 9

    for i in eachindex(param_names)
        @info param_names[i]

        ##### I loglikelihood calculation
        median_val = param_vals[i]
        param_range = median_val .* [0.1, 0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
        param_mat = hcat(fill(param_vals, nvals)...)
        param_mat[i, :] = param_range
        @time ll = ll_batch(param_mat; seed)

        ##### II plotting
        nplots = size(ll, 2)
        mean_logps = vec(mean(ll; dims = 2))

        fig = Figure(; resolution = (500, 800))

        Axis(fig[1, 1];
            title = param_names[i],
            ylabel = "Log likelihood")
        for p in 1:nplots
            scatterlines!(param_range, ll[:, p]; color = (:grey, 0.5))
        end
        scatterlines!(param_range, mean_logps; markersize = 20, color = :red)
        scatter!([median_val], [mean_logps[nvals ÷ 2 + 1]];
            color = :blue, markersize = 25)

        Axis(fig[2, 1];
            ylabel = "Log likelihood",
            xlabel = "Parameter value")
        scatterlines!(param_range, mean_logps; markersize = 15, color = :red)
        scatter!([median_val], [mean_logps[nvals ÷ 2 + 1]];
            color = :blue, markersize = 25)

        display(fig)
        save("../tmp/$(param_names[i])_noseed.png", fig)
    end
end
