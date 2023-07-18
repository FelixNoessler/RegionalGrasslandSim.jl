module PET

using Turing
using LazyArrays
using LinearAlgebra
using JLD2
using Unitful

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
    y ~ arraydist(LazyArray(@~ truncated.(Normal.(μ, σ); lower=0)))

    return (; μ=μ, σ=σ)  # extract with generated_quantities
end


function create_cyclic_effect(d; freqs = [1,2,10], period=365)
    return [sinpi.(2 .* freqs' .* d ./ period) cospi.(2 .* freqs' .* d ./ period)]
end



function predict_pet(; nyears=1, explo, datapath)
    posterior_chain = load(
        "$datapath/scenarios/evapotranspiration.jld2",
        explo)

    t_pred = repeat(1:365, nyears)
    cyclic_pred = create_cyclic_effect(t_pred; freqs = [1,2,3])
    m_pred = evaporation_model(cyclic_pred)

    pred_dat = Array(
        Turing.predict(
            m_pred,
            sample(posterior_chain, 1))
    )

    return vec(pred_dat' .* u"mm / d")
end

end # of module
