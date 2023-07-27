module Traits

using Unitful
using Distributions
using DataFrames
using JLD2
using LinearAlgebra

struct GM
    μ::Any
    Σ::Any
    ϕ::Any
end

include("traits_output.jl")

function load_data(datapath)
    ########### parameters for gaussian mixture model
    μ, Σ, ϕ = load("$datapath/input/traits_gaussian_mixture.jld2", "μ", "Σ", "ϕ")
    global gm = GM(μ, Σ, ϕ)

    return nothing
end

function inverse_logit(x)
    return exp(x) / (1 + exp(x))
end

function random_traits(n; back_transform = true)
    m = MixtureModel([
            MvNormal(gm.μ[1, :], Hermitian(gm.Σ[1, :, :])),
            MvNormal(gm.μ[2, :], Hermitian(gm.Σ[2, :, :]))],
        gm.ϕ)

    log_logit_traits = rand(m, n)

    if back_transform
        transformations = [
            exp, exp, exp, exp, exp,
            inverse_logit,
            exp, exp, exp,
        ]

        traits = Array{Float64}(undef,
            n,
            length(transformations))
        traitdf_names = [
            "LA_log", "LFM_log", "LDM_log",
            "BA_log", "SRSA_log", "AMC_logit",
            "CH_log",
            "LDMPM_log", "LNCM_log",
        ]
        trait_names = first.(split.(traitdf_names, "_"))

        for (i, t) in enumerate(transformations)
            trait = t.(log_logit_traits[i, :])
            traits[:, i] .= trait
        end

        trait_df = DataFrame(traits, trait_names)
        trait_df.SLA = trait_df.LA ./ trait_df.LDM
        trait_df.SLA = trait_df.SLA ./ 1000

        trait_df.SRSA_above = trait_df.SRSA .* trait_df.BA

        return trait_df
    else
        return log_logit_traits
    end
end

function relative_traits(; trait_data)
    trait_data = ustrip.(trait_data)
    ntraits = size(trait_data, 2)

    #### calculate extrema from more data
    many_traits = random_traits(100;)
    many_traits = Matrix(ustrip.(many_traits))

    for i in 1:ntraits
        mint, maxt = quantile(many_traits[:, i], [0.05, 0.95])
        trait_data[:, i] .= (trait_data[:, i] .- mint) ./ (maxt - mint)
    end

    [trait_data[trait_data[:, i] .< 0, i] .= 0.0 for i in 1:ntraits]
    [trait_data[trait_data[:, i] .> 1, i] .= 1.0 for i in 1:ntraits]

    return trait_data
end

end
