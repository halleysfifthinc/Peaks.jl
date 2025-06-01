```@setup benchmark_plots
include("plots/benchmark_plots.jl")
```

# Peak finding benchmarks

Peaks.jl functions can be very fast, but some features (wider window widths, support for
`missing`, etc) come at a performance penalty. Below, 3 example peak finding functions[^1]
and [`Images.findlocalmaxima`](https://juliaimages.org/latest/function_reference/#Images.findlocalmaxima)
are benchmarked against [`argmaxima`](@ref) and [`simplemaxima`](@ref) from Peaks.jl. (Only
maxima type functions are benchmarked below, but similar results can be found for
equivalent minima finding functions.)

[^1]: Currently/previously used by other Julia packages or public code. Those packages/codes
    are intentionally not named because they are not primarily peak finding packages, and
    the functions aren't part of their public API.

One important difference between the Peaks.jl functions (`argmaxima`, `simplemaxima`, etc)
and `Images.findlocalmaxima` or the example functions is that only Peaks.jl functions
recognize plateaus. Depending on the type and origins of your data recognizing plateaus may
or may not be an important/relevant feature. (Quantized data, such as data sampled from
physical measurements using ADC's, has a higher likelihood of plateaus.) As a result,
Peaks.jl functions are doing more work, but are still competitive with other obvious
implementations.

```@example benchmark_plots
times_plot # hide
```

The data generation for the above benchmark was a sine function with a peak every 5
elements; similar results are observed for more sparsely spaced peaks. However, that
benchmark may not accurately reflect real-world performance when peaks are spaced more
variably[^2].

[^2]: This benchmark will be distorted by how well the
    [CPU branch-predictor](https://en.wikipedia.org/wiki/Branch_predictor) can learn each
    function; more predictable branches will make the function run faster than could be
    expected with more realistic data.

When benchmarking against randomly sampled data, we see roughly the same performance
characteristics, now presented relative to the speed of `simplemaxima` (higher is slower):

```@example benchmark_plots
ratio_plot # hide
```

The large difference in speed between `simplemaxima` and other functions is due to the use
of SIMD code for suitable input arrays and element types; see the docstring for its
limitations compared to `argmaxima`.

**Benchmarked functions**:

!!! details "`naivemaxima`"
    ```julia
    function naivemaxima(x)
       pks = Int[]
       i = firstindex(x) + 1
       for i in firstindex(x)+1:lastindex(x)-1
           if x[i+1] < x[i] > x[i-1]
               push!(pks, i)
           end
       end
       return pks
    end
    ```

!!! details "`diffmaxima`"
    ```julia
    function diffmaxima(x)
        dx = diff(x)
        return findall((diff(dx) .< 0) .& (dx[begin:end-1] .> 0) .& (dx[begin+1:end] .< 0)) .+ firstindex(x)
    end
    ```

!!! details "`viewsmaxima`"
    ```julia
    function viewsmaxima(x)
        x1 = @view x[begin:end - 2]
        x2 = @view x[begin+1:end - 1]
        x3 = @view x[begin+2:end]
        return axes(x,1)[2:end-1][x1 .< x2 .> x3]
    end
    ```

