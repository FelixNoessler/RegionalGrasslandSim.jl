# Growth

```@meta
CurrentModule = RegionalGrasslandSim.Growth
```

the net growth of the plants is modelled by...
- the [potential growth](@ref pot_growth) that is multiplied by some [growth reducer functions](@ref reducer_functions) and a [belowground competition function](@ref below_competition), these processes are included in the main function [`growth`](@ref)
- [Leaf senescence](@ref)
- [Agricultural defoliation](@ref)

```@docs
growth
```

---
## [Potential growth](@id pot_growth)

```@docs
potential_growth
calculate_LAI
```

----
## [Reducer functions](@id reducer_functions)
The growth of each plant species in each patch is dependent on... 
- ☀ the photosynthetically active radiation [`radiation_reduction`](@ref)
- the height of the plants in relation to the community weighted mean height [`height_influence`](@ref)
- 🌡 the air temperature [`temperature_reduction`](@ref)
- 💧 the [soil water content](@ref water_stress)
- the [plant-available nutrients](@ref nut_stress)
- 📈 a seasonal effect, that is modelled by the accumulated degree days [`seasonal_reduction`](@ref)


```@docs
radiation_reduction
height_influence
community_weighted_mean_height
temperature_reduction
water_reduction
nutrient_reduction
seasonal_reduction
```
--

## [Below-ground competition](@id below_competition)

```@docs
below_ground_competition
similarity_matrix
```
--- 
## Leaf senescence

```@docs
senescence
seasonal_component_senescence
```

---
## Agricultural defoliation

Biomass is removed by...
- 🐄 [`grazing`](@ref) and [`trampling`](@ref)
- 🚜 [`mowing`](@ref)


```@docs
grazing
mowing
trampling
```
--- 
