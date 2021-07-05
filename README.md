# Peaks.jl

[![version](https://juliahub.com/docs/Peaks/version.svg)](https://juliahub.com/ui/Packages/Peaks/3TWUM)
[![pkgeval](https://juliahub.com/docs/Peaks/pkgeval.svg)](https://juliahub.com/ui/Packages/Peaks/3TWUM)
[![CI](https://github.com/halleysfifthinc/Peaks.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/halleysfifthinc/Peaks.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/halleysfifthinc/Peaks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/halleysfifthinc/Peaks.jl)
![Maintenance](https://img.shields.io/maintenance/yes/2021)



Peaks.jl contains peak (local extrema) finding functions for vector data. Contributions welcome.

## Functions

- `argmaxima`/`argminima`
  - Find the indices of the local extrema of `x` where each extrema is
    either the maximum of `x[-w:w]` or the first index of a plateau.
    If `strictbounds` is `true`, all elements of `x[-w:w]` must exist
    and may not be `missing` or `NaN`. If `strictbounds` is `false`,
    elements of `x[-w:w]` may not exist (eg peaks may be less than `w`
    indices from either end of `x`), or may be `missing` or `NaN`.
  - Supports OffsetArrays
  - See docstring for more information

- `findmaxima`/`findminima` => (indices, values)
  - Return the indices and values of local extrema

- `peakprom`
  - Find all local extrema and peak prominences in `x` matching the
    conditions `w` and `minprom`. `w` sets the minimum allowed distance
    between extrema. `minprom` sets the minimum prominence (inclusive) of
    returned extrema.
    Peak prominence is calculated as the difference between the current
    extrema and the most extreme of the smallest extrema of the lower and upper
    bounds. Bounds extend from the current extrema to the next element
    more extreme than the current extrema, or the end of the signal,
    which ever comes first.
  - See docstring for more information

## Related

- [**Images.jl**](https://github.com/JuliaImages/Images.jl)
  - [`findlocalmaxima`](https://juliaimages.org/stable/function_reference/#Images.findlocalmaxima)/[`findlocalminima`](https://juliaimages.org/stable/function_reference/#Images.findlocalminima)
    - Supports more than 1 dimension
    - Doesn't support `missing`, different window sizes
