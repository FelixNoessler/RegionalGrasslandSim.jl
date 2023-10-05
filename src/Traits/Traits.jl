module Traits

using Unitful
using Distributions
using JLD2
using LinearAlgebra
import Random

function load_gm(datapath)
    ########### parameters for gaussian mixture model
    μ, Σ, ϕ = load("$datapath/input/traits_gaussian_mixture.jld2", "μ", "Σ", "ϕ")

    m = MixtureModel([
            MvNormal(μ[1, :], Hermitian(Σ[1, :, :])),
            MvNormal(μ[2, :], Hermitian(Σ[2, :, :]))],
        ϕ)

    global m

    return nothing
end

function inverse_logit(x)
    return exp(x) / (1 + exp(x))
end

function random_traits(n; seed)
    Random.seed!(seed)

    log_logit_traits = rand(m, n)

    transformations = [
        exp, exp, exp, exp, exp,
        inverse_logit,
        exp, exp, exp,
    ]

    units = [
        u"mm^2", u"mg", u"mg",      # leaf traits
        NoUnits, u"m^2/g", NoUnits, # root traits,
        u"m",                       # LEDA
        u"g/g", u"mg/g",            # TRY
    ]

    traits = []
    traitdf_names = [
        "LA_log", "LFM_log", "LDM_log",
        "BA_log", "SRSA_log", "AMC_logit",
        "height_log",
        "LDMPM_log", "LNCM_log",
    ]
    trait_names = Symbol.(first.(split.(traitdf_names, "_")))

    for (i, t) in enumerate(transformations)
        trait = t.(log_logit_traits[i, :])
        unit_vector = repeat([units[i]], n)
        push!(traits, trait .* unit_vector)
    end

    trait_dict = Dict(trait_names .=> traits)
    trait_dict[:SLA] = trait_dict[:LA] ./ trait_dict[:LDM]
    trait_dict[:SLA] = uconvert.(u"m^2/g", trait_dict[:SLA])
    trait_dict[:SRSA_above] = trait_dict[:SRSA] .* trait_dict[:BA]

    return trait_dict
end

function relative_traits(; trait_data, large_trait_data)
    rel_traits_dict = Dict()

    for trait in keys(trait_data)
        mint, maxt = quantile(large_trait_data[trait], [0.05, 0.95])

        rel_trait_vals = (trait_data[trait] .- mint) ./ (maxt - mint)
        rel_trait_vals[rel_trait_vals .< 0.0] .= 0.0
        rel_trait_vals[rel_trait_vals .> 1.0] .= 1.0

        rel_traits_dict[trait] = rel_trait_vals
    end

    return rel_traits_dict
end

function generate_random_traits(n; seed)
    trait_dict = random_traits(n; seed)
    many_trait_dict = random_traits(100; seed)

    rel_traits = relative_traits(;
        trait_data = trait_dict, large_trait_data = many_trait_dict)

    return (; trait_dict...), (; rel_traits...)
end

@doc raw"""
    similarity_matrix(; scaled_traits, similarity_exponent)

Computes the trait similarity of all plant species.

The trait similarity between plant species $i$ and
plant species $u$ for $T$ traits is calculated as follows:
```math
\text{trait_similarity}_{i,u} =
    1-\frac{\sum_{t=1}^{t=T}
        |\text{scaled_trait}_{t,i} - \text{scaled_trait}_{t,u}|}{T}
```

To give each functional trait an equal influence,
the trait values have been scaled by the 5 % ($Q_{0.05, t}$)
and 95 % quantile ($Q_{0.95, t}$) of trait values of 100 plant species:
```math
\text{scaled_trait}_{t,i} =
    \frac{\text{trait}_{t,i} - Q_{0.05, t}}
    {Q_{0.95, t} - Q_{0.05, t}}
```

If the rescaled trait values were below zero or above one, the values were
set to zero or one respectively.
"""
function similarity_matrix(; scaled_traits, belowtrait_similarity_exponent)
    nspecies, ntraits = size(scaled_traits)
    similarity_mat = Array{Float64}(undef, nspecies, nspecies)

    for i in 1:nspecies
        diff = zeros(nspecies)

        for n in 1:ntraits
            species_trait = scaled_traits[i, n]
            diff .+= abs.(species_trait .- scaled_traits[:, n])
        end

        similarity_mat[i, :] .= 1 .- (diff ./ ntraits)
    end

    ### include here the competition strength
    return similarity_mat .^ belowtrait_similarity_exponent .* u"ha / kg"
end

end
