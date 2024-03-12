"""
    peakwidths(indices, x, proms; [strict=true, relheight=0.5, min, max]) -> (indices, widths, ledge, redge)
    peakwidths(pks::NamedTuple; [strict=true, relheight=0.5, min, max]) -> NamedTuple

Calculate the widths of peak `indices` in `x` at a reference level based on `proms` and
`relheight`, and removing peaks with widths less than `min` and/or greater than
`max`. Returns the peaks, widths, and the left and right edges at the reference level.

Peak width is the distance between the signal crossing a reference level before and after
the peak. Signal crossings are linearly interpolated between indices. The reference level is
the difference between the peak height and `relheight` times the peak prominence. Width
cannot be calculated for a `NaN` or `missing` prominence.

If a NamedTuple `pks` is given, a new NamedTuple is returned with filtered copies of fields
from `pks`. `pks` must have `:indices`, `:heights`, and `:proms` fields. If `pks` has
`:widths` and `:edges` fields, they will not be recalculated, but filtered only. Any
remaining fields will be copied unmodified.

If `strict == true`, the width for a peak with a gap in the signal (e.g. `NaN`, `missing`)
at the reference level will match the value/type of the signal gap. Otherwise, the signal
crossing will be linearly interpolated between the edges of the gap.

See also: [`peakwidths!`](@ref), [`peakproms`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = Float64[0,5,2,2,3,3,1,4,0];

julia> pks = findmaxima(x) |> peakproms!(;max=2);

julia> peakwidths(pks)
(indices = [5], heights = [3.0], data = [0.0, 5.0, 2.0, 2.0, 3.0, 3.0, 1.0, 4.0, 0.0], proms = [1.0], widths = [1.75], edges = [(4.5, 6.25)])

julia> x[4] = NaN;

julia> peakwidths(pks.indices, x, pks.proms)
([5], [NaN], [NaN], [6.25])

julia> peakwidths(pks.indices, x, pks.proms; strict=false)
([5], [2.25], [4.0], [6.25])
```
"""
function peakwidths(
    peaks::AbstractVector{Int}, x::AbstractVector, proms::AbstractVector;
    strict=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
    min=minwidth, max=maxwidth
)
    if !isnothing(minwidth)
        Base.depwarn("Keyword `minwidth` has been renamed to `min`", :peakwidths)
    end
    if !isnothing(maxwidth)
        Base.depwarn("Keyword `maxwidth` has been renamed to `max`", :peakwidths)
    end
    if !isnothing(min) || !isnothing(max)
        _peaks = copy(peaks)
    else
        # peaks will not be modified
        _peaks = peaks
    end
    peakwidths!(_peaks, x, proms; strict=strict, relheight=relheight,
        min=min, max=max)
end

peakwidths(pks::NamedTuple; kwargs...) = peakwidths!(deepcopy(pks); kwargs...)

"""
    peakwidths(; [strict, relheight, min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that calculates and filters the peak widths of a
copy of its argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0,5,2,3,3,1,4,0]) |> peakproms() |> peakwidths(; min=1.5)
(indices = [4], heights = [3], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[1], widths = Union{Missing, Float64}[1.75], edges = Tuple{Union{Missing, Float64}, Union{Missing, Float64}}[(3.5, 5.25)])
```
"""
peakwidths(; kwargs...) = function _curried_peakwidths(pks)
    return peakwidths(pks; kwargs...)
end

function _inner_widthscalcloop!(op::O, cmp::C, x::AbstractVector{T}, peaks::AbstractVector{Int}, proms::AbstractVector{Union{Missing,T}}, ledge::AbstractVector{V}, redge::AbstractVector{V}, relheight::U, _bad, fst, lst, strict::Bool) where {O,C,T,V,U}
    for i in eachindex(peaks, ledge, redge)
        prom = proms[i]
        if ismissing(prom) || isnan(prom)
            redge[i] = _bad
            ledge[i] = _bad
        else
            ht = op(x[peaks[i]], relheight * prom)
            lo = findprev(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])
            hi = findnext(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])

            if !strict
                if !isnothing(lo)
                    lo1 = findnext(v -> !ismissing(v) && cmp(ht, v), x, lo + 1)
                    @assert !isnothing(lo1)
                    lo += (ht - x[lo]) / (x[lo1] - x[lo]) * (lo1 - lo)
                end
                if !isnothing(hi)
                    hi1 = findprev(v -> !ismissing(v) && cmp(ht, v), x, hi - 1)
                    @assert !isnothing(hi1)
                    hi -= (ht - x[hi]) / (x[hi1] - x[hi]) * (hi - hi1)
                end
            else
                !isnothing(lo) && (lo += (ht - x[lo]) / (x[lo+1] - x[lo]))
                !isnothing(hi) && (hi -= (ht - x[hi]) / (x[hi-1] - x[hi]))
            end
            redge[i] = something(hi, lst)
            ledge[i] = something(lo, fst)
        end
    end

    return nothing
end

"""
    peakwidths!(indices, x; [strict=true, relheight=0.5, min, max]) -> (indices, widths, ledge, redge)
    peakwidths!(pks::NamedTuple; [strict=true, relheight=0.5, min, max]) -> NamedTuple

Calculate the widths of peak `indices` in `x` at a reference level based on `proms` and
`relheight`, removing peaks with widths less than `min` and/or greater than `max`.
Returns the modified peaks, widths, and the left and right edges at the reference level.

If a NamedTuple `pks` is given, a new NamedTuple is returned with the same fields (references)
from `pks`. `pks` must have `:indices`, `:heights`, and `:proms` fields. If `pks` has
`:widths` and `:edges` fields, they will not be recalculated, but filtered only. Any
remaining fields will be copied unmodified.

See also: [`peakwidths`](@ref), [`peakproms`](@ref), [`findmaxima`](@ref)
#
# Examples
```jldoctest ; filter = r"(\\d*)\\.(\\d{3})\\d*" => s"\\1.\\2***"
julia> x = Float64[0,5,2,2,3,3,1,4,0];

julia> pks = findmaxima(x) |> peakproms!();

julia> peakwidths!(pks; min=1)
(indices = [2, 5], heights = [5.0, 3.0], data = [0.0, 5.0, 2.0, 2.0, 3.0, 3.0, 1.0, 4.0, 0.0], proms = [5.0, 1.0], widths = [1.333, 1.75], edges = [(1.5, 2.833), (4.5, 6.25)])

julia> peakwidths!(pks.indices, pks.data, pks.proms; min=1)
([2, 5], [1.333, 1.75], [1.5, 4.5], [2.833, 6.25])
```
"""
function peakwidths!(
    peaks::AbstractVector{Int}, x::AbstractVector{T}, proms::AbstractVector{U};
    strict=true, relheight=0.5, minwidth=nothing, maxwidth=nothing,
    min=minwidth, max=maxwidth
) where {T,U}
    if !isnothing(minwidth)
        Base.depwarn("Keyword `minwidth` has been renamed to `min`", :peakwidths!)
    end
    if !isnothing(maxwidth)
        Base.depwarn("Keyword `maxwidth` has been renamed to `max`", :peakwidths!)
    end
    if !isnothing(min) && !isnothing(max)
        min < max || throw(ArgumentError("Keyword `min` must be less than `max`"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))

    # if peaks was calculated with strict=false, first(peaks) could be minima at firstindex
    if ismaxima(first(peaks), x; strict=false)
        maxima = true
    elseif isminima(first(peaks), x; strict=false)
        maxima = false
    else
        throw(ArgumentError("The first peak in `indices` is not a local extrema"))
    end

    V1 = promote_type(T, U)
    _bad = Missing <: V1 ? missing : float(Int)(NaN)

    V = promote_type(V1, float(Int))
    ledge = similar(proms, V)
    redge = similar(proms, V)
    widths = similar(proms, V)

    if strict
        lst, fst = _bad, _bad
    else
        lst = lastindex(x)
        fst = firstindex(x)
    end

    if maxima
        _inner_widthscalcloop!(-, ≤, x, peaks, proms, ledge, redge, relheight, _bad, fst, lst, strict)
    else
        _inner_widthscalcloop!(+, ≥, x, peaks, proms, ledge, redge, relheight, _bad, fst, lst, strict)
    end

    widths .= redge .- ledge

    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(widths)))
        hi = something(max, typemax(Base.nonmissingtype(eltype(widths))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ hi), widths)
        deleteat!(peaks, matched)
        deleteat!(ledge, matched)
        deleteat!(redge, matched)
        deleteat!(widths, matched)
    end

    return peaks, widths, ledge, redge
end

function peakwidths!(pks::NamedTuple; strict=true, relheight=0.5, min=nothing, max=nothing)
    !haskey(pks, :proms) && throw(ArgumentError(
        "Argument `pks` is expected to have prominences (`:proms`) already calculated"))
    if xor(hasproperty(pks, :widths), hasproperty(pks, :edges))
        throw(ArgumentError("Argument `pks` is expected have neither or both of the fields `:widths` and `:edges`."))
    end
    if !hasproperty(pks, :widths)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, widths, leftedges, rightedges = peakwidths(pks.indices, pks.data, pks.proms; relheight, strict)
        pks = merge(pks, (; widths, edges=collect(zip(leftedges, rightedges))))
    end
    filterpeaks!(pks, :widths; min, max)
    return pks
end

"""
    peakwidths!(; [strict, relheight, min, max]) -> Function

Create a function, `f(pks::NamedTuple)`, that calculates and filters (mutates) the peak
widths and other fields of its argument, `pks`, using any given keyword arguments.

# Examples
```jldoctest
julia> findmaxima([0,5,2,3,3,1,4,0]) |> peakproms!() |> peakwidths!(; min=1.5)
(indices = [4], heights = [3], data = [0, 5, 2, 3, 3, 1, 4, 0], proms = Union{Missing, Int64}[1], widths = Union{Missing, Float64}[1.75], edges = Tuple{Union{Missing, Float64}, Union{Missing, Float64}}[(3.5, 5.25)])
```
"""
peakwidths!(; kwargs...) = function _curried_peakwidths!(pks)
    return peakwidths!(pks; kwargs...)
end

