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

## Finding peaks

```@setup tutorial
using Peaks, Plots; gr()
Plots.reset_defaults()

a = 3
b = 2
c = 1

fs = 100
T = 1/fs
f1 = 5
f2 = 10
f3 = 30

t = T:T:25

sinf(t) = a*sin(2*pi*f1*T*t) + b*sin(2*pi*f2*T*t) + c*sin(2*pi*f3*T*t)

x = sinf.(t);
```

```@example tutorial
p = plot(t, x; label="signal") # hide
```

To find the peaks in your data you can use the `findmaxima` function:

```@repl tutorial
indices, heights = findmaxima(x)
```

When the peaks are plotted over the data, we see that all the local maxima have been identified.

```@example tutorial
plot!(p, t[indices], heights; seriestype=:scatter, label="maxima") # hide
```

## Peak characteristics

Two commonly desired peak characteristics can be determined using the `peakproms` and `peakwidths` functions:

```@repl tutorial
indices, proms = peakproms(indices, x)

indices, widths, edges... = peakwidths(indices, x, proms)
```

## Plotting

The peaks, prominences, and widths can be visualized all together using `plotpeaks`:

```@example tutorial
plotpeaks(t, x; peaks=indices, prominences=true, widths=true)
```

