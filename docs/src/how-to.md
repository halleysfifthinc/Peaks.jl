```@setup spacing
using Peaks, Plots, DataFrames; gr()

using Random
Random.seed!(0xFEED)
```
Signals often contain peaks that you aren't interested in. This guide will show you how to
filter peaks you don't want.

## How to filter peaks by peak spacing

Real data typically has noise that can create lots of closely spaced, unwanted peaks.

```@example spacing
T = 1/25
t = 0:T:23

multisin(t) = 3sinpi(0.1t) + 2sinpi(0.2t) + sinpi(0.6t)

y = multisin.(t) .+ 0.1rand(length(t))

pks = findmaxima(y)
plotpeaks(t, pks)
```

The simplest way to remove those peaks (assuming the signal is already filtered) is by
setting the window `w` argument in `findmaxima` and friends:

```@example spacing
pks = findmaxima(y, 15)
f = plotpeaks(t, pks)
```

If only the peaks circled in blue are wanted, then setting the window `w` too wide won't
work, since there are larger peaks that would become dominant.

```@example spacing
pks = findmaxima(y, 15) # hide
wpks = peakproms(pks; max=1) # hide
plt = plot(t[wpks.indices], wpks.heights; seriestype=:scatter, markershape=:circle, label="", # hide
    markersize=10, markercolor=RGBA(1,1,1,0), markerstrokealpha=.1, markerstrokecolor=:blue, # hide
    markerstrokewidth=2, z_order=1) # hide
plotpeaks!(plt, t, pks) # hide
```

## How to filter peaks by peak characteristics

Every peak-characteristic finding function can optionally filter the newly calculated
characteristics using the keyword arguments `min` and `max`.

Plotting all the peak characteristics and/or looking at the characteristic values can help
show which characteristics should be filtered to remove all the unwanted peaks.

```@example spacing
pks = findmaxima(y, 15) # hide
wpks = peakproms(pks; max=1) # hide
plt = plot(t[wpks.indices], wpks.heights; seriestype=:scatter, markershape=:circle, label="", # hide
    markersize=10, markercolor=RGBA(1,1,1,0), markerstrokealpha=.1, markerstrokecolor=:blue, # hide
    markerstrokewidth=2, z_order=1) # hide
plotpeaks!(plt, t, pks; show_prominences=true, show_widths=true)
```
```@example spacing
pks = peakproms(pks) |> peakwidths # hide
DataFrame(pks[Not(:data)])
```

Looking at the figure and the characteristic values, we can list the usefulness of each
characteristic for filtering:

  - Peak height?
    - There are other peaks around the same height as the peaks we want, so applying a
      `min` or `max` height filter would remove peaks we want, or allow peaks we don't
      want.
  - Peak prominence?
    - All the peaks we want have similarly small prominences (<1) and the other peaks
      have much larger prominences (>2). This would be a good filtering option, using
      `peakproms(pks; max=1)`.
  - Peak width?
    - The peaks we want have fairly similar widths (~15-25 elements wide), and the other
      peaks have larger widths (>40 elements wide). This would be a good filter, using
      `peakwidths(pks; max=30)`.

In this case, filtering by peak prominence would be the better choice, because calculating
peak widths depends on prominences, so filtering by peak prominence would do the job while
avoiding unnecessary work.

In many cases, the desired peaks aren't very different from many other peaks in any one peak
characteristic. In these situations, it may be necessary to filter multiple times based on
different peak characteristics or different `min`/`max` thresholds. There is also a
[`filterpeaks!`](@ref) function which allows you to give a filter predicate and filter by multiple
characteristics at once.



