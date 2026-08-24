t=0:1/100:1
y = 2*sin.(5*t)+3*sin.(10*t)+2*sin.(30*t)
pks, vals = findmaxima(y)
_, proms = peakprom(pks, y; strict=false)
promlinesx = vec(vcat(reshape(repeat(t[pks]; inner=2),2,length(pks)), fill(NaN, (1,length(pks)))))
promlinesy = vec([vals - proms vals fill(NaN, (length(pks),1))]');

_, widths, lower, upper = peakwidth(pks, y, proms; strict=true, relheight=prevfloat(1.))
fullwidthlinesx = vec([(lower .- 1)./fs (upper .- 1)./fs fill(NaN, (length(pks),1))]')
fullwidthlinesy = vec(vcat(reshape(repeat(vals - proms.*1; inner=2),2,length(pks)), fill(NaN, (1,length(pks)))))

_, widths, lower, upper = peakwidth(pks, y, proms; strict=true)
halfwidthlinesx = vec([(lower .- 1)./fs (upper .- 1)./fs fill(NaN, (length(pks),1))]')
halfwidthlinesy = vec(vcat(reshape(repeat(vals - proms.*0.5; inner=2),2,length(pks)), fill(NaN, (1,length(pks)))))

p = plot([
    scatter(;x=t, y, line_color="black", name="Signal", line_shape="spline", legendrank=1),
    scatter(;x=promlinesx, y=promlinesy, mode="lines", line_color="#4063d8", name="Prominence", legendrank=3),
    scatter(;x=fullwidthlinesx, y=fullwidthlinesy, mode="lines", line_color="#4063d8", showlegend=false),
    scatter(;x=t[pks], y=vals, marker=attr(;color="#cb3c33", size=8), mode="markers", name="Maxima", legendrank=2),
    scatter(;x=halfwidthlinesx, y=halfwidthlinesy, line_dash="dot", line_color="grey", name="Width", legendrank=4),
    ], Layout(;
    width=600,
    height=400,
    margin=attr(;
        t=5,
        r=5,
        b=25,
        l=30,
        pad=5,
        autoexpand=false,
    ),
    font_size="16",

))
savefig(p, "../docs/src/assets/images/peaks_prom_width.png"; width=600, height=400, scale=4)
