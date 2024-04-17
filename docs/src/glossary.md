```@setup prominence
    using Peaks, PlotlyJS

    a = 3
    b = 2
    c = 1

    fs = 100
    T = 1/fs
    f1 = .05
    f2 = .10
    f3 = .30

    sinf(t) = a*sin(2*pi*f1*t) + b*sin(2*pi*f2*t) + c*sin(2*pi*f3*t)
    t2 = 10:T:23.5
    x = sinf.(t2)
    pks, = findmaxima(x)
    currpk = pks[2]
    prevpk = pks[1]
    nextpk = pks[4]
    leftmin = findnextminima(x, prevpk)
    rightmin = findnextminima(x, currpk)

    _, proms = peakproms(pks, x)
    pkprom = proms[2]

    p = plot(scatter(;x=t2, y=x), Layout(font_size=14); config=PlotConfig(responsive=true,scrollZoom=false))
    add_hline!(p, x[currpk]; line_dash="dash", line_width=2.5, line_color="rgba(89,105,112,0.40)")
    add_trace!(p, scatter(;x=t2[pks[[1,4]]], y=x[pks[[1,4]]], mode=:markers))
    add_trace!(p, scatter(;x=[t2[currpk]], y=[x[currpk]], mode=:markers, marker_color=:purple))
    add_trace!(p, scatter(;x=[t2[pks[3]]], y=[x[pks[3]]], mode=:markers, marker_color=:gray))
    add_trace!(p, scatter(;x=t2[[leftmin,rightmin]], y=x[[leftmin,rightmin]], mode=:markers, marker_color=:blue))
    relayout!(p; annotations=[
        attr(;x=22.2, y=x[currpk]+.6, ay=x[currpk], ax=0,
            ayref="y",
            yref="y",
            yanchor="top",
            showarrow=true,
            arrowhead=2,
            arrowwidth=2.5,
            arrowsize=.7,
            arrowcolor="rgba(89,105,112,0.40)",
            standoff=0,
            startstandoff=0,
            text="Minimum reference maxima",
        ),
        attr(;x=t2[prevpk], y=x[prevpk],
            text="Previous reference peak",
            ax=15,
            showarrow=true,
            arrowhead=2,
            arrowsize=1.3,
            xshift=2,
            yshift=4,
        ),
        attr(;x=t2[currpk], y=x[currpk],
            text="Current peak",
            ax=15,
            showarrow=true,
            arrowhead=2,
            arrowsize=1.3,
            xshift=2,
            yshift=4,
        ),
        attr(;x=t2[pks[3]], y=x[pks[3]],
            text="Ignored (too small)",
            hovertext="Peak is smaller than the current peak",
            showarrow=true,
            arrowhead=2,
            arrowsize=1.3,
            xshift=-2,
            yshift=4,
        ),
        attr(;x=t2[nextpk], y=x[nextpk],
            text="Next reference peak",
            ax=-15,
            showarrow=true,
            arrowhead=2,
            arrowsize=1.3,
            xshift=-2,
            yshift=4,
        ),
        attr(;x=t2[leftmin], y=x[leftmin],
            ax=-15, ay=15,
            yanchor="top",
            text="Larger minimum",
            showarrow=true,
            arrowhead=2,
            arrowsize=1.3,
            xshift=-2,
            yshift=-4,
            ),
        attr(;x=t2[currpk], y=x[leftmin], ax=t2[currpk], ay=x[currpk]-.07,
            axref="x", ayref="y",
            yanchor="bottom",
            showarrow=true,
            arrowside="end+start",
            arrowhead=2,
            arrowsize=1.3,
        ),
        attr(;x=14.4, y=x[leftmin] + (x[currpk]-.07-x[leftmin])/2,
            xanchor="left",
            yanchor="middle",
            text="Prominence: $(round(pkprom, digits=2))",
            showarrow=false,
        ),
    ],
    margin=attr(autoexpand=false, b=10, l=10, r=10, t=10),
    showlegend=false,
    )
    add_shape!(p, line(t2[leftmin], t2[currpk]+.15, x[leftmin], x[leftmin]; line_width=1.5, line_color=:blue))
```

```@setup plateau
    using Peaks, PlotlyJS

    t = 0:.001:1
    y=clamp.(sinpi.(t), 0, 0.75)

    p = plot(scatter(;x=t, y), Layout(;
        margin=attr(autoexpand=false, b=10, l=10, r=10, t=10),
        showlegend=false,
    ))
    pks = findmaxima(y)

    add_trace!(p, scatter(;x=t[pks[1]], y=y[pks[1]], mode=:markers))
```

```@setup width
    using Peaks, PlotlyJS

    t = -1.2:.001:1.2
    y = cospi.(t) .+ (t./-10)
    pks = findmaxima(y)
    pks = peakproms!(pks)
    relheight = 0.7
    pks = peakwidths!(pks; relheight)

    xneg1 = findfirst(x -> x ≈ -1, t)

    p = plot(scatter(;y), Layout(;
        margin=attr(b=10, l=10, r=10, t=10),
        xlabel="Indices",
        showlegend=false,
    ))

    add_hline!(p, pks.heights[1] - pks.proms[1]; line_dash="dash", line_width=2, line_color="rgba(89,105,112,0.40)")
    add_hline!(p, pks.heights[1]; line_dash="dash", line_width=2, line_color="rgba(89,105,112,0.40)")

    add_trace!(p, scatter(;x=pks.indices, y=pks.heights, mode=:markers))

    codefont = attr(;
            font_family="JuliaMono,SFMono-Regular,Menlo,Consolas,Liberation Mono,DejaVu Sans Mono,monospace",
            bgcolor="#f2f2f2",
            borderpad=4,
        )
    relayout!(p, font_size=13, annotations=[
        attr(;
            x=pks.edges[1][1], ax=pks.edges[1][2],
            y=pks.heights[1] - pks.proms[1]*relheight,
            ay=pks.heights[1] - pks.proms[1]*relheight,
            xanchor="left",
            axref="x", ayref="y",
            showarrow=true,
            arrowside="end+start",
            arrowhead=2,
            arrowwidth=2,
            ),
        attr(;
            x=xneg1, ax=xneg1, y=pks.heights[1] - pks.proms[1], ay=pks.heights[1],
            axref="x", ayref="y",
            yanchor="bottom",
            showarrow=true,
            arrowside="end+start",
            arrowhead=2,
            arrowwidth=2,
        ),
        attr(;
            x=xneg1, y=1,
            xanchor="left",
            yanchor="top",
            xshift=8, yshift=-6,
            text="Prominence",
            showarrow=false,
        ),
        attr(codefont;
            x=findfirst(x -> x ≈ 0.55, t),
            xanchor="left",
            y=pks.heights[1],
            yshift=-5,
            yanchor="top",
            text="relheight = 0",
            showarrow=false,
        ),
        attr(codefont;
            x=pks.edges[1][2],
            xshift=-20,
            xanchor="right",
            y=pks.heights[1] - pks.proms[1]*relheight,
            yshift=2,
            yanchor="bottom",
            text="relheight = $(relheight)",
            showarrow=false,
        ),
        attr(codefont;
# x=0.75,
            x=findfirst(x -> x ≈ 0.75, t),
            xshift=-20,
            xanchor="right",
            y=pks.heights[1] - pks.proms[1],
            yshift=5,
            yanchor="bottom",
            text="relheight = 1",
            showarrow=false,
        ),
        attr(codefont;
# x=-0.48,
            x=findfirst(x -> x ≈ -0.55, t),
            xanchor="left",
            y=pks.heights[1] - pks.proms[1]*relheight,
            yanchor="top",
            yshift=-6,
            text="width = $(round(pks.widths[1], digits=2))",
            showarrow=false,
        ),
    ],
        margin=attr(b=10, l=10, r=10, t=10),
    )
```

# Common terminology

##### Peak (a.k.a. [local] extrema, maxima, minima, etc.)

An element `x[i]` which is more extreme than its adjacent elements, or more extreme than all
elements in `x[i-w:i+w]` where `w` is a positive integer

!!! note
    May specifically refer to the index (i.e. location) of a peak, which is most
    broadly relevant when speaking of a specific peak

##### Plateau

A "flat" peak, where the value of the maxima occurs multiple times consecutively, but surrounding elements are less than the maximal value. The first maximal value is considered the peak location. Uncommon for floating-point data.

!!! note "Example plateau"
    ```@example plateau
    p # hide
    ```

##### Peak prominence

For maxima, peak prominence is the absolute difference in height between a maxima and the
larger of the two minimums in the adjacent reference ranges. Reference ranges are the range between the current maxima and adjacent (i.e. previous and next) reference maxima, which must be at least as large as the current maxima. The same is true of minima with opposite extrema/extremum (e.g. minima for maxima, and maximum for minimum, etc.).

!!! note "Example prominence calculation"
    ```@example prominence
    p # hide
    ```

##### Peak width

Peak width is measured as the distance (in units of indices) between the intersection of a horizontal reference line with the signal on either side of a peak, where the height of the reference line is offset from the peak height by a proportion of the peak prominence.

!!! note "Example peak width calculation"
    ```@example width
    p # hide
    ```

##### "Strict"-ness

The algorithms for finding peaks and calculating peak characteristics (e.g. prominence, etc.) have specific requirements of input data to return the requested outputs. All functions are `strict == true` by default (keyword argument), but setting the `strict` keyword to `false` allows participating functions to relax some alogithmic requirements. When `strict == false`, functions will make optimistic assumptions in an attempt to return useful information (e.g. not `NaN` or `missing`) when data violates normal algorithmic requirements.

!!! warning "List of Peaks.jl function behavior/assumptions for `strict == false`"
    - `maxima`/`minima` finding functions (e.g. `findmaxima`, etc.) assume that any missing data in a window is consistent with a peak. For example:
        - The maximal/minimal value in an incomplete window (i.e. the value `i` is within `w` elements of the array end, `i+w > lastindex(x)`, or beginning, `i-w < firstindex(x)`) is assumed to be a peak (i.e. if the array continued, the data would be less/more the current maximal/minimal value)
        - The maximal/minimal value in a window containing `missing` or `NaN` elements is assumed to be a peak (i.e. the `missing` or `NaN` values would be less/more than the current value if they existed or were real numbers)
    - `peakproms` uses the larger present (i.e. not `NaN` or `missing`) value of the minimum values in each reference range (see [prominence definition](#Peak-prominence))
    - `peakwidths` linearly interpolates across a gap at the width reference level

