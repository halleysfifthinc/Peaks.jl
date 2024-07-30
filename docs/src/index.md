# Peaks.jl

Peaks.jl is a library for finding [peaks](glossary.md#peak) (i.e. local maxima and minima)
and peak characteristics (e.g. how prominent or wide are peaks, etc) in vector (1D) data
signals.

## Installation

Peaks.jl can be installed from the Julia REPL by running

```julia-repl
] add Peaks
```

```@example
mktempdir() do tmpdir # hide
    out = read(`julia --color=yes -e " # hide
        using Pkg; # hide
        redirect_stderr(devnull) # hide
        Pkg.activate(\"$tmpdir\"); # hide
        redirect_stderr(stdout) # hide
        Pkg.add(\"Peaks\");"`, String); # hide
    out = replace(out, tmpdir => "~/.julia/environments/v$(VERSION.major).$(VERSION.minor)") # hide
    print(out); # hide
end; # hide
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

To find the peaks in your data you can use the [`findmaxima`](@ref) function, which returns
a [`NamedTuple`](@extref Base Core.NamedTuple) with fields for each calculated peak
characteristic (indices, heights) and a reference to the data.

```@repl tutorial
indices, heights = findmaxima(y)
```

When the peaks are plotted over the data, we see that all the local maxima have been
identified.

```@example tutorial
plot!(p, t[indices], heights; seriestype=:scatter, label="maxima") # hide
```

### Peak characteristics

Peaks have various characteristics, including [height](glossary.md#height),
[prominence](glossary.md#prominence), and [width](glossary.md#width). Peaks.jl exports
functions for finding each characteristic: [`peakheights`](@ref), [`peakproms`](@ref),
[`peakwidths`](@ref):

```@repl tutorial
indices, proms = peakproms(indices, y)

indices, widths, edges... = peakwidths(indices, y, proms)
```

Mutating bang (`'!'`) functions are available (i.e. [`peakproms!`](@ref), etc.) when
allocations are a concern.

### Peaks `NamedTuple` & pipable API

There are Peaks.jl functions that bundle the peaks, peak characteristics, and signal into a convenient `NamedTuple`:

```@repl tutorial
pks = findmaxima(y);
pks = peakproms(pks);
pks = peakwidths(pks)
```

Mutating functions are also available for the `NamedTuple` functions; the vectors within the
`NamedTuple` are mutated and re-used in the returned tuple. The `NamedTuple` functions can also be [chained/piped](@extref Base.:|>):

```@repl tutorial
pks = findmaxima(y) |> peakproms!(;strict=false) |> peakwidths!(; max=100)
```

!!! warning "Performance tip"
    Be aware that the `NamedTuple` functions allocate more memory than the functions with
    direct/explicit arguments. If maximum performance is needed, mutating functions (e.g.
    [`peakproms!`](@ref), etc) and/or the direct, non-`NamedTuple` methods are a better choice.

### Plotting

The peaks, prominences, and widths can be visualized all together using a `Plots.jl`
recipe `plotpeaks`:

```@example tutorial
using Plots
plotpeaks(t, y; peaks=indices, prominences=true, widths=true)
```

