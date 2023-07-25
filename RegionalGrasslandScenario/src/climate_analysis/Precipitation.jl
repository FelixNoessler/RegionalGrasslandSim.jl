module Precipitation

using Statistics
using StatsBase
using Random
using Distributions
using Turing
using LinearAlgebra
using LazyArrays
using JLD2
using Unitful

function load_data(datapath)
    hai_samples = load("$datapath/scenarios/prec.jld2", "HAI")
    alb_samples = load("$datapath/scenarios/prec.jld2", "ALB")
    sch_samples = load("$datapath/scenarios/prec.jld2", "SCH")

    global posterior = (;
        HAI=hai_samples,
        ALB=alb_samples,
        SCH=sch_samples,
    )

    return nothing
end

function transform_t(; t, original_t=t)
    t_min, t_max = extrema(original_t)
    return (t .- t_min) ./ (t_max - t_min)
end

##############################
struct ZeroLogNormal{T<:Real} <: ContinuousUnivariateDistribution
    ϕ₀::T
    μ::T
    σ::T
    ZeroLogNormal{T}(ϕ₀, μ, σ) where {T} = new{T}(ϕ₀, μ, σ)
end

function ZeroLogNormal(ϕ₀::T, μ::T, σ::T; check_args::Bool=true) where {T <: Real}
    Distributions.@check_args ZeroLogNormal (ϕ₀, one(ϕ₀) ≥ ϕ₀ ≥ zero(ϕ₀)) (σ, σ ≥ zero(σ))
    return ZeroLogNormal{T}(ϕ₀, μ, σ)
end

ZeroLogNormal(ϕ₀::Real, μ::Real, σ::Real; check_args::Bool=true) = ZeroLogNormal(promote(ϕ₀, μ, σ)...; check_args=check_args)
ZeroLogNormal(ϕ₀::Integer, μ::Integer, σ::Integer; check_args::Bool=true) = ZeroLogNormal(float(ϕ₀), float(μ), float(σ); check_args=check_args)

Distributions.@distr_support ZeroLogNormal 0.0 Inf

# 0 - helper function
StatsBase.params(d::ZeroLogNormal) = (d.ϕ₀, d.μ, d.σ)

# 1 - rand
function Base.rand(rng::AbstractRNG, d::ZeroLogNormal)
    ϕ₀, μ, σ = params(d)

    if rand(rng) < ϕ₀
        return 0.0
    else
        return rand(rng, truncated(LogNormal(μ, σ); upper=190))
    end
end

# 2 sampler
Distributions.sampler(rng::AbstractRNG, d::ZeroLogNormal) = Base.rand(rng::AbstractRNG, d::ZeroLogNormal)

# 3 - logpdf
function Distributions.logpdf(d::ZeroLogNormal, x::Real)
    ϕ₀, μ, σ = params(d)
    l = truncated(LogNormal(μ, σ); upper=190)

    if x == 0
        return log(ϕ₀)
    elseif x > 0
        return log(1 - ϕ₀) + logpdf(l, x)
    else
        return -Inf
    end
end

function create_cyclic_effect(d; freqs = [1,2,10], period=365)
    return [sinpi.(2 .* freqs' .* d ./ period) cospi.(2 .* freqs' .* d ./ period)]
end


invlogit(x::Real) = exp(x)/(1+exp(x))

@model function precip_model(c)
    ### priors
    αμ ~ Normal(0.0, 2)
    βμ ~ MvNormal(zeros(size(c, 2)), 0.1 * I)
    αϕ₀ ~ Normal(0.0, 0.5) # on logit scale
    βϕ₀ ~ MvNormal(zeros(size(c, 2)), 0.1 * I) # on logit scale
    σ ~ truncated(Normal(0, 3); lower=0)

    ### cyclic effects
    μ = αμ .+ c * βμ
    ϕ₀ = invlogit.(αϕ₀ .+ c * βϕ₀)

    ### likelihood
    # y ~ arraydist(LazyArray(@~ ZeroLogNormal.(ϕ₀, μ, σ)))
    y ~ arraydist(ZeroLogNormal.(ϕ₀, μ, σ))

    return (ϕ₀=ϕ₀, μ=μ)  # extract with generated_quantities
end

function predict_precipitation(;
    nyears,
    explo)

    t_pred = 1:365
    cyclic_pred = create_cyclic_effect(t_pred; freqs = [1,2,3])
    m_pred = precip_model(cyclic_pred);

    pred_matrix = Array(Turing.predict(
        m_pred,
        sample(posterior[Symbol(explo)], nyears)))
    pred_vector = reshape(pred_matrix', :)

    return pred_vector .* u"mm / d"
end

end
