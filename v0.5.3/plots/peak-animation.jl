using Peaks, PlotlyJS

include(joinpath(@__DIR__, "standards.jl"))

T = .1
t = round.(0.:T:23; digits=1)
ax = axes(t, 1)

sinf(t) = 3sinpi(0.1t) + 2sinpi(0.2t) + sinpi(0.6t)
y = round.(sinf.(t); sigdigits=3)

frame_len = 600*T
clamp_window(i) = clamp(i-w,ax):clamp(i+w,ax)

w = 16
pks = findmaxima(y, w)

pks_trace = scatter(;x=t[pks.indices], y=pks.heights, mode="markers", zorder=10,
name="Maxima", legendrank=2, marker=attr(;color="#d62728", opacity=zeros(Int, size(pks.indices))))
maximum_trace = scatter(;x=[t[argmax(y[clamp_window(1)])]], y=[maximum(y[clamp_window(1)])],
    mode="markers", marker_color="blueviolet", legendrank=3, zorder=5, name="Window maximum")
windowshade = scatter(;x=[0], y=[0], visible="legendonly", opacity=.65,
    line=attr(;color=:black, width=15), name="Window")

p = Plot([
        scatter(;x=t, y, name="Signal", legendrank=1),
        maximum_trace,
        pks_trace,
        windowshade],
    Layout(;
        font_size=14,
        legend=attr(;
            itemclick=false,
            itemdoubleclick=false,
        ),
        margin=attr(autoexpand=false, b=10, l=10, r=10, t=10),
        yaxis_fixedrange=true,
        xaxis_fixedrange=true,
        xaxis_range=[first(t),last(t)],
        shapes = [
            rect(t[clamp(1-w, ax)], t[clamp(1+w,ax)], 0, 1;
                       yref="y domain", fillcolor=:black, opacity=0.25, layer="below",
                       line_width=0),
            vline(t[clamp(1, ax)])
        ],
        updatemenus = [attr(;
            x=0, y=0,
            active=0,
            type="buttons",
            yanchor="bottom",
            xanchor="left",
            direction="left",
            pad=attr(;l=6, b=6),
            buttons=[
                attr(;
                    method="animate",
                    label="Play",
                    args=[nothing, attr(;
                        mode="next",
                        fromcurrent=true,
                        transition=attr(;
                            duration=frame_len,
                            easing="linear",
                        ),
                        frame = attr(;
                            duration=frame_len,
                            redraw=false
                        )
                )]),
                attr(;
                    method="animate",
                    label="Pause",
                    args=[[nothing], attr(;
                        mode="next",
                        fromcurrent=true,
                        transition_duration=0,
                        frame = attr(;
                            duration=0,
                            redraw=false
                        )
                )]),
                attr(;
                    method="animate",
                    label="Reset",
                    args=[[1], attr(;
                        mode="next",
                        transition_duration=0,
                        frame = attr(;
                            duration=0,
                            redraw=false
                        )
                )])
            ]
        )],
        sliders = [ attr(;
            x=1, xanchor="right",
            pad = attr(;l=-8, r=-8, t=20),
            len = 1,
            transition=attr(;
                duration=frame_len,
                easing="linear",
            ),
            currentvalue_visible=false,
            ticklen=0,
            minorticklen=0,
            steps = [
                attr(;
                method="animate",
                label="",
                args=[[i], attr(;
                    mode="immediate",
                    transition=attr(;
                        duration=frame_len,
                        easing="linear",
                    ),
                    frame = attr(;
                        redraw=false,
                        duration=0,
                    )
                )])
                for i in ax[begin:5:end]
            ]
        )]),
    [
        frame(; name = i,
            data = [
                attr(;
                    x=[t[clamp_window(i)[argmax(y[clamp_window(i)])]]],
                    y=[y[clamp_window(i)[argmax(y[clamp_window(i)])]]]
                ),
                attr(;marker_opacity=map(x -> x ≥ i ? 0 : 1, pks.indices)),
            ],
            traces = [1,2],
            layout = attr(;
                shapes = [
                    rect(t[clamp(i-w, ax)], t[clamp(i+w,ax)], 0, 1)
                    vline(t[clamp(i, ax)])
                ]))
        for i in ax
    ]; config=PlotConfig(;displayModeBar=false,showLink=false,showEditInChartStudio=false))
