## Finding peaks

```@docs
findpeaks
findmaxima
findminima
argmaxima
argminima
simplemaxima
simpleminima
maxima
minima
```

## Peak characteristics & filtering

```@docs
peakproms
peakproms!
peakwidths
peakwidths!
peakheights
peakheights!
filterpeaks!
```

## Convenience functions

```@docs
findnextmaxima
findnextminima
ismaxima
isminima
isplateau
```

## Plotting

Plotting functions are provided as package extensions. [`plotpeaks`](@ref)/[`plotpeaks!`](@ref)
require [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to be loaded (e.g. `using Plots`),
and [`peaksplot`](@ref)/[`peaksplot!`](@ref) require a
[Makie.jl](https://github.com/MakieOrg/Makie.jl) backend to be loaded (e.g. `using CairoMakie`).

```@docs
plotpeaks
plotpeaks!
peaksplot
peaksplot!
```
