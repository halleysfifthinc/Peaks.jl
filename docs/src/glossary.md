```@setup peak-animation
include("plots/peak-animation.jl")
```

```@setup prominence
include("plots/prominence-animation.jl")
```

```@setup plateau
include("plots/plateau.jl")
```

```@setup width
include("plots/width.jl")
```

# Common terminology

##### [Peak (a.k.a. [local] extrema, maxima, minima, etc.)](@id peak)

An element `x[i]` which is more extreme than its adjacent elements, or more extreme than all
elements in the window `x[i-w:i+w]` where `w` is a positive integer.

In the animation below, the maximum in the window `x[i-w:i+w]` is shown as the purple dot.
When the location of the window maximum matches the current index (the vertical black line),
a peak is identified (red dot). Use the scroll bar at the bottom to understand why some
"peaks" aren't found. (Hint: pay attention the window size and the window maximum.)

```@example peak-animation
PlotForceHTML(p) # hide
```

!!! note
    "Peak" may specifically refer to the index (i.e. location) of the peak, which is most
    broadly relevant when speaking of a specific peak

##### Plateau

A "flat" peak, where the value of the extrema occurs multiple times consecutively, but
surrounding elements are less than the extremum. The first occurence of the extrema is
considered the peak location. Uncommon for
[floating-point data](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Floating-Point-Numbers).

!!! note "Example plateau"
    ```@example plateau
    p # hide
    ```

##### [Peak prominence](@id prominence)

For maxima, peak prominence is the absolute difference in height between a maxima and the
larger of the two minimums in the adjacent reference ranges. Reference ranges cover the data
between the current maxima and adjacent (i.e. previous and next) reference maxima, which
must be at least as large as the current maxima, or the beginning/end of the array. The same
is true of minima with opposite extrema/extremum (e.g. minima for maxima, and maximum for
minimum, etc.).

```@example prominence
PlotForceHTML(p) # hide
```

##### [Peak width](@id width)

Peak width is measured as the distance (in units of indices) between the intersection of a
horizontal reference line with the signal on either side of a peak, where the height of the
reference line is offset from the peak height by a proportion of the peak prominence
(keyword argument `relheight` for the `peakwidths` functions).

!!! note "Example peak width calculation"
    ```@example width
    p # hide
    ```

##### ["Strict"-ness](@id strict)

The default behavior of peak finding and related functions (e.g. `peakproms`, etc.) is to
only return results that are exactly correct, and to return nothing (i.e. ignore a potential
peak),
[`NaN`](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Special-floating-point-values),
or [`missing`](https://docs.julialang.org/en/v1/manual/missing/), as appropriate for a given
function. This behavior is controlled by the `strict` keyword argument (`true` by default).
Setting the `strict` keyword to `false` allows these functions to relax some data
requirements. When `strict == false`, functions will make optimistic assumptions in an
attempt to return useful information (e.g. not `NaN` or `missing`) when data violates
default requirements. **This can produce results that are not technically correct, but
sometimes this is desired/needed.**

**`strict`-ness should only affect new peaks/characteristics (i.e. only peaks
detected with `strict == false`). Any observed behavior otherwise (i.e. non-`strict` peak
characteristics are altered) is undesired[^1] and an
[issue](https://github.com/halleysfifthinc/Peaks.jl/issues/new/choose) should be opened.**

[^1]: We try to ensure that non-`strict` peak characteristics will be
    unaffected during `strict == false` mode. Any issues with examples of altered
    characteristics are appreciated, and fixes will be attempted, but we do **not** guarantee
    that altered peak characteristics can be prevented.

!!! warning "List of Peaks.jl function behavior/assumptions for `strict == false`"
    - `maxima`/`minima` finding functions (e.g. `findmaxima`, etc.) assume that any missing
      data in a window is consistent with a peak. For example:
        - The maximal/minimal value in an incomplete window (e.g. an index `i` within `w`
          elements of the array beginning or end, `i-w < firstindex(x)` or `i+w >
          lastindex(x)`) is assumed to be a peak (i.e. if the array continued, the data
          would be less/more the current maximal/minimal value). This allows the first or
          last elements of an array to be considered peaks.
        - The maximal/minimal value in a window containing `missing` or `NaN` elements is
          assumed to be a peak (i.e. the `missing` or `NaN` values would be less/more than
          the current value if they existed or were real numbers)
    - `peakproms` uses the larger present (i.e. not `NaN` or `missing`) value of the minimum
      values in each reference range (see [prominence definition](#Peak-prominence))
    - `peakwidths` linearly interpolates across a gap at the width reference level

