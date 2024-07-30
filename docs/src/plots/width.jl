using Peaks, PlotlyJS

t = round.(-1.2:.01:1.2; digits=2)
y = round.(cospi.(t) .+ (t./-10); digits=3)
pks = findmaxima(y)
pks = peakproms!(pks)
relheight = 0.7
pks = peakwidths!(pks; relheight)

xneg1 = findfirst(x -> x ≈ -1, t)

p = Plot(scatter(;y), Layout(;
    margin=attr(b=10, l=10, r=10, t=10),
    xlabel="Indices",
    showlegend=false,
))

add_hline!(p, pks.heights[1] - pks.proms[1]; line_dash="dash", line_width=2, line_color="rgba(89,105,112,0.40)")
add_hline!(p, pks.heights[1]; line_dash="dash", line_width=2, line_color="rgba(89,105,112,0.40)")

add_trace!(p, scatter(;x=pks.indices.-1, y=pks.heights, mode=:markers))

codefont = attr(;
        font_family="JuliaMono,SFMono-Regular,Menlo,Consolas,Liberation Mono,DejaVu Sans Mono,monospace",
        bgcolor="#f2f2f2",
        borderpad=4,
    )
relayout!(p, font_size=13, annotations=[
    attr(;
        x=pks.edges[1][1]-1, ax=pks.edges[1][2]-1,
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
        x=pks.edges[1][2]-1,
        xshift=-20,
        xanchor="right",
        y=pks.heights[1] - pks.proms[1]*relheight,
        yshift=2,
        yanchor="bottom",
        text="relheight = $(relheight)",
        showarrow=false,
    ),
    attr(codefont;
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
