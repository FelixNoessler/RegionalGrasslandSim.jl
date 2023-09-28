module Landuse

using Dates
using Unitful

function fertilisation_input(;
    dates,
    nutrient_inputs,
    nyears)
    fertilisation = zeros(365)

    for i in eachindex(dates)
        my_date = Date(dates[i], dateformat"m-d")
        fertilisation[dayofyear(my_date)] = nutrient_inputs[i]
    end

    return repeat(fertilisation, nyears)
end

function mowing_input(;
    mowing_doys,
    nyears,
    cutting_height)
    mowing_height = zeros(365) .* u"m"

    for mowing_doy in mowing_doys
        mowing_height[mowing_doy] = cutting_height
    end

    return repeat(mowing_height, nyears)
end

function grazing_input(;
    grazing_start,
    grazing_end,
    grazing_intensity,
    nyears)
    grazing_density = fill(0.0u"ha ^ -1", 365)

    for i in eachindex(grazing_start)
        start_date = Date(grazing_start[i], dateformat"mm-dd")
        end_date = Date(grazing_end[i], dateformat"mm-dd")
        livestock_density = grazing_intensity[i]

        grazing_density[dayofyear(start_date):dayofyear(end_date)] .= livestock_density .*
                                                                      u"ha ^ -1"
    end

    return repeat(grazing_density, nyears)
end

end
