module PET

using Turing
using LinearAlgebra
using JLD2
using Unitful

function load_data(datapath)
    hai_samples = load("$datapath/scenarios/evapotranspiration.jld2", "HAI")
    alb_samples = load("$datapath/scenarios/evapotranspiration.jld2", "ALB")
    sch_samples = load("$datapath/scenarios/evapotranspiration.jld2", "SCH")

    global posterior = (;
        HAI = hai_samples,
        ALB = alb_samples,
        SCH = sch_samples)

    return nothing
end

@model function evaporation_model(c)
    ### priors
    αμ ~ Normal(0.0, 2)
    βμ ~ MvNormal(zeros(size(c, 2)), 0.1 * I)
    ασ ~ Normal(0.0, 2)
    βσ ~ MvNormal(zeros(size(c, 2)), 1 * I)

    ### cyclic effects
    μ = αμ .+ c * βμ
    σ = exp.(ασ .+ c * βσ)

    ### likelihood
    y ~ arraydist(truncated.(Normal.(μ, σ); lower = 0))

    return (; μ = μ, σ = σ)  # extract with generated_quantities
end

function create_cyclic_effect(d; freqs = [1, 2, 10], period = 365)
    return [sinpi.(2 .* freqs' .* d ./ period) cospi.(2 .* freqs' .* d ./ period)]
end

function predict_pet(; nyears = 1, explo)
    t_pred = repeat(1:365, nyears)
    cyclic_pred = create_cyclic_effect(t_pred; freqs = [1, 2, 3])
    m_pred = evaporation_model(cyclic_pred)

    pred_dat = Array(Turing.predict(m_pred,
        sample(posterior[Symbol(explo)], 1)))

    return vec(pred_dat' .* u"mm / d")
end

end # of module
