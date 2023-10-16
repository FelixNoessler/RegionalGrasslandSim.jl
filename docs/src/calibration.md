# Calibration data from the Biodiversity Exploratories project

### Input data


### Calibration data

#### Raw data
| Variable                | Description                                                                           | time span      |
| ----------------------- | ------------------------------------------------------------------------------------- | -------------- |
| `biomass[plot, year]`   | Dried aboveground biomass, cutted at a height of 4 cm once per year in spring [g m⁻²] | 2009 - 2022    |
| `soilmoisture[plot, t]` | Daily soil moisture, data availability dependent on plots [%]                         | 2009 - ongoing |


**Data sources:**\
biomass: [explo16209v2](@cite) [explo12706v2](@cite) [explo14346v3](@cite) [explo15588v2](@cite) [explo16826v4](@cite) [explo19807v4](@cite) [explo19809v3](@cite) [explo21187v3](@cite) [explo23486v4](@cite) [explo24166v4](@cite) [explo26151v4](@cite) [explo27426v5](@cite) [explo31180v22](@cite) [explo31387v10](@cite)\
soil moisture: [explo19007v6](@cite)

#### Community weighted mean traits
| Variable                      | Description                                      |
| ----------------------------- | ------------------------------------------------ |
| `CWM_sla[year, patch]`        | Specific leaf area [m² g⁻¹]                      |
| `CWM_amc[year, patch]`        | Arbuscular mycorrhizal colonisation [-]          |
| `CWM_srsa_above[year, patch]` | Root surface area / aboveground biomass [m² g⁻¹] |
| `CWM_height[year, patch]`     | Plant height [m]                                 |
| `CWM_lncm[year, patch]`       | Leaf nitrogen / leaf mass [mg g⁻¹]               |

**Data sources:**\
leaf trait data (sla): [explo24807v2](@cite)\
root trait data (srsa_above, amc): [explo26546v2](@cite)\
height and lncm data: [trydb](@cite)