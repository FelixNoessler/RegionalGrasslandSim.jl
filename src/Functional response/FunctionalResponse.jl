module FunctionalResponse

"""
    amc_nut_response(;
        mycorrhizal_colon,
        max_right_upper_bound=1,
        min_right_upper_bound=0.7,
        max_AMC_half_response=0.6,
        min_AMC_half_response=0.05,
        mid_AMC=0.2,
        slope=10
    )

Transforms the mycorrhizal colonisation into parameters
of the response curve of growth in relation to nutrient
availability.
"""
function amc_nut_response(;
    mycorrhizal_colon,
    max_right_upper_bound = 1,
    min_right_upper_bound = 0.7,
    max_AMC_half_response = 0.6,
    min_AMC_half_response = 0.05,
    mid_AMC = 0.35,
    slope = 10)

    #### check parameter
    in_range = 0.0 .<= mycorrhizal_colon .<= 1.0
    if !all(in_range)
        error("$mycorrhizal_colon (mycorrhizal_colonisation) not between 0 and 1")
    end

    #### Denominator for two logistic functions that
    #### transforms the mycorrizal colonisation to parameters
    #### of the reponse curve of the growth to water availability
    denominator = @. (1 + exp(-slope * (mycorrhizal_colon - mid_AMC)))

    #### right upper bound of repsonse curve
    bounds_diff = (min_right_upper_bound - max_right_upper_bound)
    response_right_bound = @. max_right_upper_bound + bounds_diff / denominator

    #### x vlaue of the midpoint of the response curve
    bounds_diff = (min_AMC_half_response - max_AMC_half_response)
    response_midpoint = @. max_AMC_half_response + bounds_diff / denominator

    return response_right_bound, response_midpoint
end

"""
    srsa_response(;
        SRSA, # specific root surface area / above ground biomass
        mid_SRSA = 0.12,
        slope_func_parameters = 20,
        min_right_upper_bound = 0.7,
        max_right_upper_bound = 1,
        min_SRSA_half_response = 0.05,
        max_SRSA_half_response = 0.6
        )

"""
function srsa_response(;
    SRSA_above, # specific root surface area / above ground biomass [m^2 g^-1]
    mid_SRSA_above = 0.12, # m^2 g^-1
    slope_func_parameters = 40,
    min_right_upper_bound = 0.7,
    max_right_upper_bound = 1,
    min_SRSA_half_response = 0.05,
    max_SRSA_half_response = 0.6)
    denominator = @. (1 + exp(-slope_func_parameters * (SRSA_above - mid_SRSA_above)))

    K = @. max_right_upper_bound +
           (min_right_upper_bound - max_right_upper_bound) / denominator
    x0 = @. max_SRSA_half_response +
            (min_SRSA_half_response - max_SRSA_half_response) / denominator

    return K, x0
end

"""
    sla_water_response(;
        SLA,
        mid_SLA = 0.025,
        slope_func_parameter = 75,
        min_SLA_half_response = -0.8,
        max_SLA_half_response = 0.8
    )

"""
function sla_water_response(;
    SLA,
    mid_SLA = 0.025,
    slope_func_parameter = 75,
    min_SLA_half_response = -0.8,
    max_SLA_half_response = 0.8)
    x0 = @. min_SLA_half_response +
            (max_SLA_half_response - min_SLA_half_response) /
            (1 + exp(-slope_func_parameter * (SLA - mid_SLA)))

    return x0
end

end
