using Peaks, PlotlyJS

t = 0:.001:1
y=clamp.(sinpi.(t), 0, 0.7)

p = make_subplots(;cols=2, shared_yaxes="rows", horizontal_spacing=0.1,
    subplot_titles=["Plateau" "Not plateaus"])
relayout!(p,
    margin=attr(autoexpand=true, b=10, l=10, r=10, t=25),
    showlegend=false,
)

pks = findmaxima(y)
add_trace!(p, scatter(;x=t, y); col=1)
add_trace!(p, scatter(;x=t[pks[1]], y=y[pks[1]], mode=:markers); col=1)

t_time = searchsortedfirst(t, 0.35):searchsortedlast(t, 0.65)
clamp_begin = t[findfirst(==(0.7), y)]
clamp_end = t[findlast(==(0.7), y)]
y[t_time] += sinpi.(range(clamp_begin, clamp_end; length=length(t_time))) .- 0.7
add_trace!(p, scatter(;x=t, y=y, line_color="#636efa"); col=2)
add_trace!(p, scatter(;x=[clamp_begin, clamp_end], y=[0.7,0.7], mode="markers",
    marker_color=:gray); col=2)
