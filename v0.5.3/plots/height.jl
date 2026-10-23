using Peaks, PlotlyJS

include(joinpath(@__DIR__, "standards.jl"))

T = .1
w = 16

t = round.(0.:T:23; digits=3)

sinf(t) = 3sinpi(0.1t) + 2sinpi(0.2t) + sinpi(0.6t)
y = round.(sinf.(t); sigdigits=3)

pks = findmaxima(y, w)

pks_trace = scatter(;x=t[pks.indices], y=pks.heights, mode="markers", zorder=10,
    yhoverformat=".2g",
    name="Maxima", legendrank=2, marker=attr(;color="#d62728"))

p = Plot([
        scatter(;x=t, y, name="Signal", legendrank=1, hoverinfo="none"),
        pks_trace
    ],
    Layout(;
        hoverdistance=100,
        font_size=14,
        legend=attr(;
            itemclick=false,
            itemdoubleclick=false,
        ),
        margin=attr(autoexpand=false, b=40, l=10, r=10, t=10),
        yaxis_fixedrange=true,
        xaxis_fixedrange=true,
        xaxis_range=[first(t),last(t)],
    ); config=PlotConfig(;displayModeBar=false,showLink=false,showEditInChartStudio=false))

