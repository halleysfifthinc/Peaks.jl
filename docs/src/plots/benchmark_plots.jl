using BenchmarkTools, PlotlyJS, DataFrames, DataFramesMeta, Statistics

include(joinpath(@__DIR__, "standards.jl"))
# include(joinpath(@__DIR__, "../../../benchmark/benchmarks.jl"))

function_colors = Dict(
    ["simplemaxima", "argmaxima", "findlocalmaxima",
     "diffmaxima", "viewsmaxima", "naivemaxima"] .=>
     PlotlyJS.PlotlyBase.templates[:plotly].layout.colorway[1:6])

results = BenchmarkTools.load(
    joinpath(@__DIR__, "../../../benchmark/vs_other_functions.json"))[1];

sparsedf = DataFrame(mapreduce(vcat, keys(results), values(results)) do f,fv
    mapreduce(vcat, filter(((k,v),) -> occursin(r"^sparse", k), fv)) do (sp, spv)
        map(pairs(spv)) do (l, v)
            spm = match(r"(?<==).*", sp).match
            (f, parse(Int, spm), parse(Int, l), v)
        end
    end
end, ["function", "sparsity", "length", "benchmark"])
sort!(sparsedf, [:length, :function])

randdf = DataFrame(mapreduce(vcat, results) do (f,fv)
    mapreduce(vcat, filter(((k,v),) -> occursin(r"^(uniform|normal)", k), fv)) do (rt, rtv)
        map(pairs(rtv)) do (l, v)
            m = match(r"^(uniform|normal)_(\d+)bit", rt)
            (f, m[1], parse(Int, m[2]), parse(Int, l), v)
        end
    end
end, ["function", "type", "bits", "length", "benchmark"])
sort!(randdf, [:length, :function])


times_df = @chain sparsedf begin
    @subset :sparsity .== 4
    @transform @byrow(:by_length = begin
                _b = mean(:benchmark)
                _b.time = (_b.time/:length)
                _b.gctime = (_b.gctime/:length)
                _b.memory = round(Int, _b.memory/:length)
                _b
            end)
    @by([:sparsity, :length],
        :function, :by_length,
        :time_ratio = BenchmarkTools.ratio.(time.(:by_length),
            Ref(time(:by_length[findfirst(==("simplemaxima"),:function)]))),
        :memory_ratio = BenchmarkTools.ratio.(memory.(:by_length),
            Ref(memory(:by_length[findfirst(==("simplemaxima"),:function)]))),
        :allocs_ratio = BenchmarkTools.ratio.(allocs.(:by_length),
            Ref(allocs(:by_length[findfirst(==("simplemaxima"),:function)]))))
end
sort!(times_df, :length);

base_layout = attr(;
    font_size=14,
    legend=attr(;
        # itemclick=false,
        # itemdoubleclick=false,
        tracegroupgap=5
    ),
    margin=attr(autoexpand=false, b=10, l=10, r=10, t=10),
    )

times_plot = Plot([
        scatter(;x=df[!, :length], y=time.(df[!, :by_length]),
            name=df[!, :function][1],
            marker_color=function_colors[df[!, :function][1]])
        for df in groupby(times_df, :function)
        ],
    Layout(;
        base_layout...,
        yaxis=attr(;
            type="log",
            title_text="Time (ns)",
        ),
        xaxis=attr(;
            type="log",
            title_text="Vector length (n)"
        )
))

ratio_df = @chain randdf begin
    @subset(:type .== "normal", :bits .== 24)
    @transform @byrow(:by_length = begin
            _b = mean(:benchmark)
            _b.time = (_b.time/:length)
            _b.gctime = (_b.gctime/:length)
            _b.memory = round(Int, _b.memory/:length)
            _b
        end)
    @by([:length, :bits, :type], :function,
        :time_ratio = BenchmarkTools.ratio.(time.(:by_length),
            Ref(time(:by_length[findfirst(==("simplemaxima"),:function)]))),
        :memory_ratio = BenchmarkTools.ratio.(memory.(:by_length),
            Ref(memory(:by_length[findfirst(==("simplemaxima"),:function)]))),
        :allocs_ratio = BenchmarkTools.ratio.(allocs.(:by_length),
            Ref(allocs(:by_length[findfirst(==("simplemaxima"),:function)]))))
end
sort!(ratio_df, :length);

trace_legend = Dict( (func, bits, trace) => if trace == :time_ratio
        bits == 8 ? true : "legendonly"
    else
        false
    end
    for func in  ["simplemaxima", "argmaxima", "findlocalmaxima",
        "diffmaxima", "viewsmaxima", "naivemaxima"]
    for bits in [8,16,20,24,32,64]
    for trace in (:time_ratio, :memory_ratio, :allocs_ratio)
)

function get_and_disable(d, k)
    v = d[k]
    if v == true
        d[k] = false
    end
    return v
end

trace_keys = [
    (df[!, :function][1], df[!, :bits][1], trace)
    for df in groupby(@subset(ratio_df, :type .== "normal"), [:bits, :function])
    for trace in (:time_ratio, :memory_ratio, :allocs_ratio)
]

ratio_traces = [
    scatter(;x=df[!, :length], y=df[!, trace],
        name=df[!, :function][1],
        marker_color=function_colors[df[!, :function][1]])
    for df in groupby(ratio_df, [:bits, :function])
    for trace in (:time_ratio,)
]

ratio_plot = Plot(ratio_traces,
    Layout(;
        base_layout...,
        yaxis=attr(;
            type="log",
            title_text="Ratio",
        ),
        xaxis=attr(;
            type="log",
            title_text="Vector length (n)"
        )
    )
)
