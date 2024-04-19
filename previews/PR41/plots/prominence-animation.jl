using Peaks, PlotlyJS

struct PlotForceHTML
    p::Plot
end

# "require-loaded" here because the first peak-animation plot already loaded plotly
Base.show(io, ::MIME"text/html", p::PlotForceHTML) = PlotlyJS.PlotlyBase.to_html(io, p.p;
    autoplay=false, include_plotlyjs="require-loaded", include_mathjax=missing,
    full_html=false);

a = 3
b = 2
c = 1

T = .1
f1 = .05
f2 = .10
f3 = .30

sinf(t) = a*sin(2*pi*f1*t) + b*sin(2*pi*f2*t) + c*sin(2*pi*f3*t)
t = 10:T:23.5
y = sinf.(t)
pks, = findmaxima(y)
currpk = pks[2]
prevpk = pks[1]
ignpk = pks[3]
nextpk = pks[4]
leftmin = findnextminima(y, prevpk)
rightmin = findnextminima(y, currpk)

_, proms = peakproms(pks, y)
pkprom = proms[2]

min_ref_max_arrow = attr(;x=21.8, y=y[currpk]+.6, ay=y[currpk], ax=0,
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
)
prev_peak_arr = attr(;x=t[prevpk], y=y[prevpk],
    text="Previous reference peak",
    ax=11.5,
    axref="x",
    showarrow=true,
    arrowhead=2,
    arrowsize=1.3,
    standoff=4,
)
curr_peak_arr = attr(;x=t[currpk], y=y[currpk],
    text="Current peak",
    ax=15, # px
    showarrow=true,
    arrowhead=2,
    arrowsize=1.3,
    standoff=4,
)
ign_peak_arr = attr(;x=t[pks[3]], y=y[pks[3]],
    text="Ignored",
    hovertext="Peak is smaller than the current peak",
    showarrow=true,
    arrowhead=2,
    arrowsize=1.3,
    standoff=4,
)
next_peak_arr = attr(;x=t[nextpk], y=y[nextpk],
    text="Next reference peak",
    ax=-15,
    showarrow=true,
    arrowhead=2,
    arrowsize=1.3,
    standoff=4,
)
ref_mins_arr = [
    attr(;x=t[leftmin], y=y[leftmin], ax=13.5, ay=-3.4,
        ayref="y", axref="x",
        yanchor="top",
        showarrow=true,
        arrowhead=2,
        arrowsize=1.3,
        standoff=5,
        text="Reference minimums",
    );
    attr(;x=t[rightmin], y=y[rightmin], ax=13.7, ay=-3.9,
        ayref="y", axref="x",
        yanchor="bottom",
        showarrow=true,
        arrowhead=2,
        arrowsize=1.3,
        standoff=12,
        startstandoff=12,
    )
]

larger_ref_min = attr(;x=t[leftmin], y=y[leftmin], ax=13.5, ay=-3.4,
        ayref="y", axref="x",
        yanchor="top",
        showarrow=true,
        arrowhead=2,
        arrowsize=1.3,
        standoff=5,
        text="Larger minimum",
)

prom_arrow = attr(;
    x=t[currpk], y=y[leftmin],
    ax=t[currpk], ay=y[currpk]-.07,
    axref="x", ayref="y",
    yanchor="bottom",
    showarrow=true,
    arrowside="end+start",
    arrowhead=2,
    arrowsize=1.3,
    text="",
)
prom_text = attr(;x=14.4, y=y[leftmin] + (y[currpk]-.07-y[leftmin])/2,
    xanchor="left",
    yanchor="middle",
    text="Prominence: $(round(pkprom, digits=2))",
    showarrow=false,
)

frame_currpeak = [
    frame(;name=0, # reset traces (color, visibility, etc)
        data = [
            attr(;
                visible=true,
                x=fill(t[currpk], 3),
                y=fill(y[currpk], 3),
                marker_color=["rgba(214,39,40,0)", "rgba(128,128,128,0)", "rgba(214,39,40,0)"],
            ),
            attr(;
                visible=true,
                x=fill(t[currpk], 2),
                y=fill(y[currpk], 2),
                marker_color=["rgba(38,93,128,0)", "rgba(38,93,128,0)"]
            ),
            attr(;
                visible=true,
                x=fill(t[leftmin], 2),
                y=fill(y[leftmin], 2),
                marker_color="rgba(38,93,128,0)",
            ),
        ],
        layout = attr(;
            annotations = nothing,
            shapes = nothing,
        ),
        traces = [1,2,3],
    ),
    frame(;name=1,
        layout = attr(;
            annotations = [
                curr_peak_arr;
            ],
            shapes = [],
        ),
    ),
    frame(;name=2,
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
            ],
            shapes = [
                hline(y[currpk]; line_dash="dash", line_width=2.5,
                    line_color="rgba(89,105,112,0.40)")
            ],
        ),
    )
]

frame_refpeaks = [
    frame(;name=3,
        data = [
            attr(;
                x=[t[prevpk]; fill(t[currpk], 2)],
                y=[y[prevpk]; fill(y[currpk], 2)],
                marker_color=["rgba(214,39,40,1)", "rgba(128,128,128,0)", "rgba(214,39,40,0)"],
            ),
        ],
        traces = [1],
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
            ],
        ),
    ),
    frame(;name=4,
        data = [
            attr(;
                x=[t[prevpk]; t[ignpk]; t[ignpk]],
                y=[y[prevpk]; y[ignpk]; y[ignpk]],
                marker_color=["rgba(214,39,40,1)", "rgba(128,128,128,1)", "rgba(214,39,40,0)"],
            ),
        ],
        traces = [1],
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
                ign_peak_arr;
            ],
        ),
    ),
    frame(;name=5,
        data = [
            attr(;
                x=[t[prevpk]; t[ignpk]; t[nextpk]],
                y=[y[prevpk]; y[ignpk]; y[nextpk]],
                marker_color=["rgba(214,39,40,1)", "rgba(128,128,128,1)", "rgba(214,39,40,1)"],
            ),
        ],
        traces = [1],
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
                ign_peak_arr;
                next_peak_arr;
            ],
        ),
    ),
]

frame_minpks = [
    frame(;name=6,
        data = [
            attr(;
                x=[t[leftmin], t[rightmin]],
                y=[y[leftmin], y[rightmin]],
                marker_color=["rgba(38,93,128,1)", "rgba(38,93,128,1)"]
            ),
        ],
        traces = [2],
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
                ign_peak_arr;
                next_peak_arr;
                ref_mins_arr;
            ],
        ),
    ),
    frame(;name=7,
        traces=[2],
        data = [
            attr(;
                x=[t[leftmin], t[rightmin]],
                y=[y[leftmin], y[rightmin]],
                marker_color=["rgba(38,93,128,1)", "rgba(38,93,128,0)"]
            ),
        ],
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
                ign_peak_arr;
                next_peak_arr;
                larger_ref_min;
                attr(;showarrow=false,text="");
            ],
        ),
    ),
]

frame_prom = [
    frame(;name=8,
        data = [
            attr(;
                marker_color="rgba(38,93,128,1)",
            ),
        ],
        traces = [3],
    ),
    frame(;name=9,
        data = [
            attr(;
                x=[t[leftmin]; 14],
            ),
        ],
        traces = [3],
    ),
    frame(;name=10,
        layout = attr(;
            annotations = [
                curr_peak_arr;
                min_ref_max_arrow;
                prev_peak_arr;
                ign_peak_arr;
                next_peak_arr;
                larger_ref_min;
                prom_text;
                prom_arrow;
            ],
        ),
    ),
]

p = Plot([
        scatter(;x=t, y, hoverinfo="x+y");
        scatter(;x=fill(t[currpk], 3), y=fill(y[currpk], 3), mode=:markers, zorder=10,
            hoverinfo="x+y", marker=attr(;
                color=["rgba(214,39,40,0)", "rgba(128,128,128,0)", "rgba(214,39,40,0)"]
            ))
        scatter(;x=fill(t[currpk], 2), y=fill(y[currpk], 2), mode=:markers, zorder=1,
            hoverinfo="x+y", marker=attr(;
                color=["rgba(38,93,128,0)", "rgba(38,93,128,0)"]
            ))
        scatter(;x=fill(t[leftmin], 2), y=fill(y[leftmin],2), marker_size=1, zorder=1000,
            marker_color="rgba(38,93,128,0)", hoverinfo="none",)
        scatter(;x=[t[currpk]], y=[y[currpk]], mode=:markers, zorder=100, marker_color=:purple,
            hoverinfo="x+y")
    ], Layout(;
        font_size=14,
        showlegend=false,
        legend=attr(;
            itemclick=false,
            itemdoubleclick=false,
        ),
        margin=attr(autoexpand=false, b=10, l=10, r=10, t=10),
        yaxis_fixedrange=true,
        xaxis_fixedrange=true,
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
                    args=[[0:10;], attr(;
                        mode="immediate",
                        transition=[
                            attr(;duration=0, easing="linear",),
                            attr(;duration=100, easing="linear",),
                            attr(;duration=100, easing="linear",),
                            attr(;duration=800, ordering="traces first"),
                            attr(;duration=800, ordering="traces first"),
                            attr(;duration=800, ordering="traces first"),
                            attr(;duration=800, ordering="traces first"),
                            attr(;duration=800, ordering="traces first"),
                            attr(;duration=200, ordering="traces first"),
                            attr(;duration=400, easing="linear",),
                            attr(;duration=800),
                        ],
                        frame = [
                            attr(;duration=10, redraw=true),
                            attr(;duration=1500, redraw=false),
                            attr(;duration=1500, redraw=false),
                            attr(;duration=2000, redraw=false),
                            attr(;duration=2000, redraw=false),
                            attr(;duration=2000, redraw=false),
                            attr(;duration=2000, redraw=false),
                            attr(;duration=1000, redraw=false),
                            attr(;duration=200, redraw=false),
                            attr(;duration=800, redraw=false),
                            attr(;duration=800, redraw=false),
                        ],
                )]),
                # attr(;
                #     label="Reset",
                #     method="update",
                #     args=[
                #         attr(;visible=false,), # data
                #         attr(; # layout
                #             annotations = nothing,
                #             shapes = nothing,
                #         ),
                #         [1,2,3] # traces
                #     ])
            ]
        )],
    ), [
        frame_currpeak;
        frame_refpeaks;
        frame_minpks;
        frame_prom;
    ]; config=PlotConfig(displayModeBar=false,showLink=false,showEditInChartStudio=false))

