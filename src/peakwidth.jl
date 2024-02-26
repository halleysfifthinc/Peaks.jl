"""
    peakwidths(peaks, x, proms;
        strict=true,
        relheight=0.5,
        minwidth=nothing,
        maxwidth=nothing
    ) -> (peaks, widths, leftedge, rightedge)

Calculate the widths of `peaks` in `x` at a reference level based on `proms` and
`relheight`, and removing peaks with widths less than `minwidth` and/or greater than
`maxwidth`. Returns the peaks, widths, and the left and right edges at the reference level.

Peak width is the distance between the signal crossing a reference level before and after
the peak. Signal crossings are linearly interpolated between indices. The reference level is
the difference between the peak height and `relheight` times the peak prominence. Width
cannot be calculated for a `NaN` or `missing` prominence.

The width for a peak with a gap in the signal (e.g. `NaN`, `missing`) at the reference level
will match the value/type of the signal gap if `strict == true`. For `strict ==
false`, the signal crossing will be linearly interpolated between the edges of the gap.

See also: [`peakproms`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)

# Examples
```jldoctest
julia> x = [0,1,0,-1.];

julia> xpks = argmaxima(x)
1-element Vector{Int64}:
 2

julia> peakwidths(xpks, x, [1])
([2], [1.0], [1.5], [2.5])

julia> x[3] = NaN;

julia> peakwidths(xpks, x, [1])
([2], [NaN], [1.5], [NaN])

julia> peakwidths(xpks, x, [1]; strict=false)
([2], [1.0], [1.5], [2.5])
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

"""
    peakwidths!(peaks, x, proms;
        strict=true,
        relheight=0.5,
        minwidth=nothing,
        maxwidth=nothing
    ) -> (peaks, widths, leftedge, rightedge)

Calculate the widths of `peaks` in `x` at a reference level based on `proms` and
`relheight`, removing peaks with widths less than `minwidth` and/or greater than `maxwidth`.
Returns the modified peaks, widths, and the left and right edges at the reference level.

See also: [`peakwidths`](@ref), [`peakproms`](@ref), [`findminima`](@ref), [`findmaxima`](@ref)
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
        min < max || throw(ArgumentError("max width must be greater than min width"))
    end
    all(∈(eachindex(x)), peaks) ||
        throw(ArgumentError("peaks contains invalid indices to x"))

    # if peaks was calculated with strict=false, first(peaks) could be minima at firstindex
    fp = length(peaks) > 1 ? peaks[2] : first(peaks)
    if fp > 1 && ((x[fp] < x[fp-1]) === true)
        pktype = :minima
    else
        pktype = :maxima
    end
    cmp = pktype === :maxima ? (≤) : (≥)
    op = pktype === :maxima ? (-) : (+)

    V1 = promote_type(T, U)
    _bad = Missing <: V1 ? missing : float(Int)(NaN)

    V = promote_type(V1, typeof(_bad))
    ledge = similar(proms, V)
    redge = similar(proms, V)

    if strict
        lst, fst = _bad, _bad
    else
        lst = lastindex(x)
        fst = firstindex(x)
    end

    for i in eachindex(peaks, ledge, redge)
        prom = proms[i]
        if ismissing(prom) || isnan(prom)
            redge[i] = _bad
            ledge[i] = _bad
        else
            ht = op(x[peaks[i]], relheight * proms[i])
            lo = findprev(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])
            up = findnext(v -> !ismissing(v) && cmp(v, ht), x, peaks[i])

            if !strict
                if !isnothing(lo)
                    lo1 = findnext(v -> !ismissing(v) && cmp(ht, v), x, lo + 1)
                    lo += (ht - x[lo]) / (x[lo1] - x[lo]) * (lo1 - lo)
                end
                if !isnothing(up)
                    up1 = findprev(v -> !ismissing(v) && cmp(ht, v), x, up - 1)
                    up -= (ht - x[up]) / (x[up1] - x[up]) * (up - up1)
                end
            else
                !isnothing(lo) && (lo += (ht - x[lo]) / (x[lo+1] - x[lo]))
                !isnothing(up) && (up -= (ht - x[up]) / (x[up-1] - x[up]))
            end
            redge[i] = something(up, lst)
            ledge[i] = something(lo, fst)
        end
    end

    widths::Vector{V} = redge - ledge

    if !isnothing(min) || !isnothing(max)
        lo = something(min, zero(eltype(widths)))
        up = something(max, typemax(Base.nonmissingtype(eltype(widths))))
        matched = findall(x -> !ismissing(x) && !(lo ≤ x ≤ up), widths)
        deleteat!(peaks, matched)
        deleteat!(ledge, matched)
        deleteat!(redge, matched)
        deleteat!(widths, matched)
    end

    return peaks, widths, ledge, redge
end

"""
    peakwidths!(pks) -> NamedTuple
    peakwidths!() -> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`.
- `max`: Filter out any peak with a height greater than `min`.
- `relheight`: How far up to measure width, as fraction of the peak prominence. Defaults to `0.5`.
- `strict`: How to handle `NaN` and `missing` values. See documentation for more details. Default to `true`.

Find the widths of the peaks in `pks`, and filter out any peak
with a width smaller than `min` or greater than `max`.
The widths are returned in the field `:widths` of the returned named tuple.
The edges of the peaks are also added in the field `:edges`.

If the positional argument `pks` is omitted, a function is returned such
that `peakwidths!(pks)` is equivalent to `pks |> peakwidths!`.

Note: If `pks` does not have a field `proms`, it is added. This is
because it is needed to calculate the peak width.

Note: This function mutates the vectors stored in the NamedTuple `pks`,
and not the named tuple itself.

See also: [`peakproms!`](@ref), [`peakheights!`](@ref)

# Examples
```jldoctest
julia> data = [1, 5, 1, 3, 2];

julia> pks = findmaxima(data);

julia> pks = peakwidths!(pks)
(indices = [2, 4], heights = [5, 3], data = [1, 5, 1, 3, 2], proms = Union{Missing, Int64}[4, 1], widths = [1.0, 0.75], edges = [(1.5, 2.5), (3.75, 4.5)])

julia> data |> findmaxima |> peakwidths!
(indices = [2, 4], heights = [5, 3], data = [1, 5, 1, 3, 2], proms = Union{Missing, Int64}[4, 1], widths = [1.0, 0.75], edges = [(1.5, 2.5), (3.75, 4.5)])
```
"""
function peakwidths!(pks::NamedTuple; minwidth=nothing, maxwidth=nothing, min=minwidth, max=maxwidth, relheight=0.5, strict=true)
    if !isnothing(minwidth)
        Base.depwarn("Keyword `minwidth` has been renamed to `min`", :peakwidths!)
    end
    if !isnothing(maxwidth)
        Base.depwarn("Keyword `maxwidth` has been renamed to `max`", :peakwidths!)
    end
    if !hasproperty(pks, :proms)  # Add proms if needed
        pks = peakproms!(pks; strict)
    end
    if xor(hasproperty(pks, :widths), hasproperty(pks, :edges))
        throw(ArgumentError("The named tuple `pks` (first argument to `peakwidths!` is expected have both the fields `:widths` and `:edges`, or to have neither of them. The provided `pks` only has one of them."))
    end
    if !hasproperty(pks, :widths)
        # Avoid filtering by min/max/strict here, so that it always happens outside if-statement.
        # Pro: one less edge case. Con: More internal allocations
        _, widths, leftedges, rightedges = peakwidths(pks.indices, pks.data, pks.proms; relheight, strict)
        pks = merge(pks, (; widths, edges=collect(zip(leftedges, rightedges))))
    end
    filterpeaks!(pks, min, max, :widths)
    return pks
end
peakwidths!(; kwargs...) = pks -> peakwidths!(pks; kwargs...)

"""
    peakwidths(pks) -> NamedTuple
    peakwidths() -> Function

# Optional keyword arguments
- `min`: Filter out any peak with a height smaller than `min`.
- `max`: Filter out any peak with a height greater than `min`.
- `relheight`: How far up to measure width, as fraction of the peak prominence. Defaults to `0.5`.
- `strict`: How to handle `NaN` and `missing` values. See documentation for more details. Default to `true`.

Non-mutation version of `peakwidths!`. Note that
this copies all vectors in `pks`, including the data.
This means that it is less performant. See the docstring for
`peakwidths!` for more information.
"""
peakwidths(pks::NamedTuple; kwargs...) = peakwidths!(deepcopy(pks); kwargs...)
