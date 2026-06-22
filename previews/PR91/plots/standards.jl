using PlotlyBase
using PlotlyBase: _TRACE_TYPES

struct PlotForceHTML
    p::Plot
end

const INCLUDE_PLOTLYJS = Ref("require-loaded")

# "require-loaded" here because the first peak-animation plot already loaded plotly
function Base.show(io, ::MIME"text/html", p::PlotForceHTML)
    PlotlyBase.to_html(io, filter_template!(p.p; no_colorscales=true);
        autoplay=false, full_html=false,
        include_plotlyjs=INCLUDE_PLOTLYJS[], include_mathjax=missing)
end

# Default template in PlotlyBase includes templates for many different trace types,
# filtering unnecessary elements can save 4.8 KiB per plot
function filter_template!(p::Plot; no_colorscales=false)
    template = deepcopy(p.layout.template)
    tracetypes = unique!(map(x -> x.type, p.data))
    rmtraces = setdiff(_TRACE_TYPES, tracetypes)
    foreach(rmtraces) do trace
        delete!(template.data, trace)
    end

    if no_colorscales
        delete!(template.layout.fields, :colorscale)
    end

    p.layout.template = template
    return p
end

