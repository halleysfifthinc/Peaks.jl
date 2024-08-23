# Peaks.jl

[![version](https://juliahub.com/docs/General/Peaks/stable/version.svg)](https://juliahub.com/ui/Packages/General/Peaks)
[![pkgeval](https://juliahub.com/docs/General/Peaks/stable/pkgeval.svg)](https://juliahub.com/ui/Packages/General/Peaks)
[![stable-docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://halleysfifthinc.github.io/Peaks.jl/stable)
[![dev-docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://halleysfifthinc.github.io/Peaks.jl/dev)
[![CI](https://github.com/halleysfifthinc/Peaks.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/halleysfifthinc/Peaks.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/halleysfifthinc/Peaks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/halleysfifthinc/Peaks.jl)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Peaks.jl contains peak (local extrema) finding functions for vector data. Visit the documentation for a complete introduction, how-to guide, and reference. Contributions welcome.

[![signal-with-peaks-prominences-and-widths](docs/src/assets/images/maxima_prom_width.png)](#)

## Features

- Find peak (maxima or minima) locations, height, prominence, and width
    - Filter peaks by peak spacing (window size), height, prominence, and width (including "Full Width Half Maximum (FWHM)")
- Fully supports `NaN`/`missing` with optional tolerance using keyword arg `strict`:
    - Conventional handling/propagation of `NaN`/`missing` when `strict = true` (the default)
    - Reasonable alternatives when `strict = false`

## Related

- [**Images.jl**](https://github.com/JuliaImages/Images.jl)
  - [`findlocalmaxima`](https://juliaimages.org/stable/function_reference/#Images.findlocalmaxima)/[`findlocalminima`](https://juliaimages.org/stable/function_reference/#Images.findlocalminima)
    - Supports more than 1 dimension
    - Doesn't support `missing`, different window sizes
