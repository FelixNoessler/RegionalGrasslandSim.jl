module AirTemperature

using JLD2
using Distributions
using Unitful


function load_data(datapath)
    hai_samples = load("$datapath/scenarios/air_temp.jld2", "HAI")
    alb_samples = load("$datapath/scenarios/air_temp.jld2", "ALB")
    sch_samples = load("$datapath/scenarios/air_temp.jld2", "SCH")
    freqs = load("$datapath/scenarios/air_temp.jld2", "freqs")

    global posterior = (;
        HAI=hai_samples,
        ALB=alb_samples,
        SCH=sch_samples,
        freqs=freqs
    )

    return nothing
end

function create_cyclic_effect(
    d;
    freqs,
    period=365)

    return [sinpi.(2 .* freqs' .* d ./ period) cospi.(2 .* freqs' .* d ./ period)]
end

function temp_ar_predict(; y_start, tend, θ, c)
    ncyclic = size(c, 2)

    αμ = θ[1]
    βμ = θ[2:ncyclic+1]
    σ = θ[ncyclic+2]
    ϕ = θ[ncyclic+3]
    ν = θ[ncyclic+4]

    ### cyclic effects
    μ = αμ .+ c * βμ
    y = [y_start]

    for t in 2:tend
        nu = (1-ϕ)*μ[t] + ϕ*y[t-1]
        y_new = rand(TDist(ν) * σ + nu)
        push!(y, y_new)
    end

    return y
end

function predict_temperature(; nyears, explo)
    tend = nyears * 365

    posterior_pred = temp_ar_predict(;
        θ=vec(Array(sample(posterior[Symbol(explo)],1))),
        tend,
        c = create_cyclic_effect(1:tend; freqs=posterior[:freqs]),
        y_start = 0.0
    )

    return posterior_pred .* u"°C"
end

"""
    yearly_temp_cumsum(d)

Calculates the temperature sum for each year.
"""
function yearly_temp_cumsum(d)
    adj = 365
    nyears = length(d) ÷ adj + 1

    final_cumsum = Array{Float64}(undef, length(d))

    d = ustrip.(d)
    d[d .< 0.0] .= 0.0

    for y in 1:nyears
        sliced_d = d[1+adj*(y-1):min(adj*y, length(d))]
        final_cumsum[1+adj*(y-1):min(adj*y, length(d))] .= cumsum(sliced_d)
    end

    return final_cumsum
end

end  # end of module
