module PAR

using Statistics
using Unitful
using Turing
using LinearAlgebra
using JLD2


function load_data(datapath)
    hai_samples = load("$datapath/scenarios/par_modeling.jld2", "HAI")
    alb_samples = load("$datapath/scenarios/par_modeling.jld2", "ALB")
    sch_samples = load("$datapath/scenarios/par_modeling.jld2", "SCH")
    hai_t = load("$datapath/scenarios/par_modeling.jld2", "HAI_t")
    alb_t = load("$datapath/scenarios/par_modeling.jld2", "ALB_t")
    sch_t = load("$datapath/scenarios/par_modeling.jld2", "SCH_t")

    global posterior = (;
        HAI=hai_samples,
        ALB=alb_samples,
        SCH=sch_samples,
        HAI_t=hai_t,
        ALB_t=alb_t,
        SCH_t=sch_t
    )

    return nothing
end

@model function par_model(c)
    βc ~ MvNormal(zeros(size(c, 2)), 20 * I)
    σ ~ truncated(Normal(0, 1.5); lower=0)
    cyclic = c * βc

    μ = cyclic
    y ~ arraydist(truncated.(Normal.(μ, σ^2), 0, 15.1))
    return (; cyclic)
end

function transform_t(; t, original_t=t)
    t_min, t_max = extrema(original_t)
    return (t .- t_min) ./ (t_max - t_min)
end

function create_cyclic_effect(d; freqs = [1], period=1 )
    return [sinpi.(2 .* freqs' .* d ./ period) cospi.(2 .* freqs' .* d ./ period)]
end

function make_prediction(; t_pred=1:365, original_t, posterior_sample)
    x_pred = transform_t(; t=t_pred, original_t)
    cyclic_predict = create_cyclic_effect(x_pred)
    y_samples = Turing.predict(par_model(cyclic_predict), posterior_sample)
    return Array(y_samples)
end

function predict_par(; explo, nyears)
    t_pred = repeat(1:365, nyears)
    posterior_sample = sample(posterior[Symbol(explo)], 1)
    Y_pred = make_prediction(;
        t_pred,
        original_t=posterior[Symbol("$(explo)_t")],
        posterior_sample)

    return vec(Y_pred) .* u"MJ / d / m^2"
end

end # of module
