function derive_WHC_PWP_nutrients(;
    patch_xdim, patch_ydim,
    nutheterog,
    sand, silt, clay, organic, bulk, rootdepth,
    total_N, CN_ratio)

    l = NeutralLandscapes.PlanarGradient(; direction = 100)
    gradientm = nothing
    if isone(patch_xdim) && isone(patch_ydim)
        gradientm = fill(0.5, 1, 1)
    else
        gradientm = rand(l, patch_xdim, patch_ydim)
    end

    WHC, PWP = input_WHC_PWP(; sand, silt, clay, organic, bulk, rootdepth)
    nutmat = input_nutrients(; gradientm, nutheterog, total_N, CN_ratio)

    WHCmat = fill(WHC, patch_xdim * patch_ydim)u"mm"
    PWPmat = fill(PWP, patch_xdim * patch_ydim)u"mm"

    return WHCmat, PWPmat, nutmat
end


@doc raw"""
    input_nutrients(; gradientm, total_N, CN_ratio)

Derive a nutrient index by combining total nitrogen and carbon to nitrogen ratio.

```math
    N = \frac{1}{1 + exp(-β₁ ⋅ \text{total_N} -β₂ ⋅ \text{CN_ratio}⁻¹)}
```

- `CN_ratio`: carbon to nitrogen ratio [-]
- `total_N`: total nitrogen [g kg⁻¹]
"""
function input_nutrients(; gradientm, nutheterog, total_N, CN_ratio)
    #### data from the biodiversity exploratories
    # mintotal_N = 1.2525
    # maxtotal_N = 30.63
    # minCN_ratio = 9.0525
    # maxCN_ratio = 13.6025

    mintotal_N = 0.0
    maxtotal_N = 50.0
    minCN_ratio = 5.0
    maxCN_ratio = 25.0

    β₁ = β₂ = 0.1
    β₃ = 10.0  * nutheterog

    N = (total_N  - mintotal_N) / (maxtotal_N - mintotal_N)
    CN = (CN_ratio - minCN_ratio) / (maxCN_ratio - minCN_ratio)

    m = @. 1 / (1 + exp(-β₁ * N -β₂ * 1/CN -β₃ * (gradientm - 0.5)))

    return vec(m)
end


@doc raw"""
    input_WHC_PWP(; sand, silt, clay, organic, bulk, rootdepth)

Derive walter holding capacity (WHC) and
permanent wilting point (PWP) from soil properties.

```math
\begin{align}
    θ₁ &= a₁ ⋅ \text{sand} + b₁ ⋅ \text{silt} + c₁ ⋅ \text{clay} +
            d₁ ⋅ \text{organic} + e₁ ⋅ \text{bulk} \\
    \text{WHC} &= θ₁ ⋅ \text{rootdepth} \\
    θ₂ &= a₂ ⋅ \text{sand} + b₂ ⋅ \text{silt} + c₂ ⋅ \text{clay} +
            d₂ ⋅ \text{organic} + e₂ ⋅ \text{bulk} \\
    \text{PWP} &= θ₂ ⋅ \text{rootdepth}
\end{align}
```

Equation and coefficients are taken from [Gupta1979](@cite).
The coefficients a, b, c, d and e differ for the water holding
capacity (matrix potential Ψ = -0.07 bar) and
the permanent wilting point (matrix potential Ψ = -15 bar).
The empirical coefficients that were estimated by [Gupta1979](@cite)
can be seen in the folling table:

| Ψ [bar] | a        | b        | c        | d        | e       |
| ------- | -------- | -------- | -------- | -------- | ------- |
| -0.07   | 0.005678 | 0.009228 | 0.009135 | 0.006103 | -0.2696 |
| -15     | -5.9e-5  | 0.001142 | 0.005766 | 0.002228 | 0.02671 |

- `sand`: sand content [%]
- `silt`: silt content [%]
- `clay`: clay content [%]
- `bulk`: bulk density [g cm⁻³]
- `organic`: organic matter content [%]
- `rootdepth`: rooting depth [mm]
- `θ`: water content [cm³ cm⁻³]
- `WHC`: water holding capacity [mm]
- `PWP`: permanent wilting point [mm]
"""
function input_WHC_PWP(; sand, silt, clay, organic, bulk, rootdepth)
    θ₁ = 0.005678 * sand + 0.009228 * silt + 0.009135 * clay +
        0.006103 * organic - 0.2696 * bulk
    WHC = θ₁ * rootdepth

    θ₂ = -5.9e-5 * sand + 0.001142 * silt + 0.005766 * clay +
        0.002228 * organic + 0.02671 * bulk
    PWP = θ₂ * rootdepth

    return WHC, PWP
end
