# Peaks.jl

## Installation

Peaks.jl can be installed from the Julia REPL by running

```julia-repl
] add Peaks
```

```@example
_, io_err = mktemp(); # hide
pkgout = redirect_stdio(stdout=devnull, stderr=io_err) do; # hide
    run(`julia --color=yes -E 'using Pkg; Pkg.activate(;temp=true); Pkg.add("Peaks");'`) # hide
end; # hide
seekstart(io_err); # hide
out = read(IOContext(io_err, :color => true), String); # hide
close(io_err); # hide
print(out); # hide
```

## Getting started

### Finding peaks

```@setup tutorial
using Peaks, Plots; gr()

T = 1/25
t = 0:T:23

f(t) = 3sinpi(0.1t) + 2sinpi(0.2t) + sinpi(0.6t)

y = f.(t);
```

```@example tutorial
p = plot(t, y; label="signal") # hide
```

To find the peaks in your data you can use the `findmaxima` function:

```@repl tutorial
indices, heights = findmaxima(y)
```

When the peaks are plotted over the data, we see that all the local maxima have been identified.

```@example tutorial
plot!(p, t[indices], heights; seriestype=:scatter, label="maxima") # hide
```

### Peak characteristics

Two commonly desired peak characteristics can be determined using the `peakproms` and `peakwidths` functions:

```@repl tutorial
indices, proms = peakproms(indices, y)

indices, widths, edges... = peakwidths(indices, y, proms)
```

Mutating bang (`'!'`) functions are available for `peakproms` (e.g. `peakproms!`),
`peakwidths`, and `peakheights`.

### Peaks `NamedTuple` & pipable API

There are Peaks.jl functions that bundle the peaks, peak characteristics, and signal into a convenient `NamedTuple`:

```@repl tutorial
pks = findmaxima(y);
pks = peakproms(pks);
pks = peakwidths(pks)
```

Mutating functions are also available for the `NamedTuple` functions; the vectors within the
`NamedTuple` are mutated and re-used in the returned tuple. The `NamedTuple` functions can also be piped:

```@repl tutorial
pks = findmaxima(y) |> peakproms!(;strict=false) |> peakwidths!(; max=100)
```

!!! warning "Performance tip"
    Be aware that the `NamedTuple` functions allocate more memory than the functions with
    direct/explicit arguments. If maximum performance is needed, mutating functions (e.g.
    [`peakproms!`](@ref), etc) and/or the direct/non-`NamedTuple` functions are a better choice.

### Plotting

The peaks, prominences, and widths can be visualized all together using the `Plots.jl`
recipe `plotpeaks`:

```@example tutorial
using Plots
plotpeaks(t, y; peaks=indices, prominences=true, widths=true)
```

