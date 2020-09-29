struct Maxima; end
struct Minima; end

for (comps, Extrema) in (((:>=, minimum, max), Maxima),
                       ((:<=, maximum, min), Minima))

    comp, extremum, extrem = comps

    @eval begin
        function peakprom(
            E::$Extrema, x::AbstractVector{T}, w=1;
            strictbounds=true, minprom=nothing
        ) where T
            if ($Extrema) === Maxima
                m = argmaxima(x, w; strictbounds=strictbounds)

                # The extremum search space in the bounding intervals can be reduced by
                # restricting the search space to known peaks/reverse peaks. The cost of
                # finding all peaks/reverse peaks should be mitigated by the fact that
                # the same peaks/reverse peaks will be the pivotal elements for
                # numerous peaks.
                if !strictbounds
                    m′ = (w === 1) ? m : argmaxima(x, 1; strictbounds=false)
                    notm = argminima(x, 1; strictbounds=false)
                end
            else
                m = argminima(x, w; strictbounds=strictbounds)
                if !strictbounds
                    m′ = (w === 1) ? m : argminima(x, 1; strictbounds=false)
                    notm = argmaxima(x, 1; strictbounds=false)
                end
            end

            M = lastindex(m)
            proms = similar(x,M)

            if strictbounds
                lbegin, lend = firstindex(x), lastindex(x)

                @inbounds for i in eachindex(m, proms)
                    # Find left and right bound (self-intersections)
                    lb = something(findprev(y -> ($comp)(y, x[m[i]]) === true, x, m[i] - 2),
                        lbegin)
                    rb = something(findnext(y -> ($comp)(y, x[m[i]]) === true, x, m[i] + 2),
                        lend)
                    # @show m[i] m[lb] m[rb]

                    # Find extremum of left and right bounds
                    if isempty(lb:(m[i] - 1))
                        lref = missing
                    else
                        lref = ($extremum)(view(x, lb:(m[i] - 1)))
                    end

                    if isempty((m[i] + 1):rb)
                        rref = missing
                    else
                        rref = ($extremum)(view(x, (m[i] + 1):rb))
                    end

                    proms[i] = abs(x[m[i]] - ($extrem)(lref, rref))
                end
            else
                m′val = x[m′]
                notmval = x[notm]

                j = something(findfirst(y -> y === m[1], m′), firstindex(m′))
                k = something(findfirst(>(m[1]), notm), firstindex(notm))
                @inbounds for i in eachindex(m, proms)
                    j = something(findnext(y -> y === m[i], m′, j), j)
                    k = something(findnext(>(m[i]), notm, k), lastindex(notm))

                    # Find left and right bounding peaks
                    _lb = something(findprev(($comp)(m′val[j]), m′val, j - 1), firstindex(m′))
                    _rb = something(findnext(($comp)(m′val[j]), m′val, j + 1), lastindex(m′))

                    # Find left and right reverse peaks just inside the bounding peaks
                    lb = something(findprev(<(m′[_lb]+2), notm, k-1), firstindex(notm))
                    rb = something(findnext(>(m′[_rb]-2), notm, k), lastindex(notm))

                    # @show m′[_lb], m′[_rb]
                    # @show notm[lb], notm[rb]

                    if isempty(lb:(k-1))
                        lref = missing
                    else
                        lref = ($extremum)(view(notmval, lb:(k - 1)))
                    end

                    if isempty(k:rb)
                        rref = missing
                    else
                        rref = ($extremum)(view(notmval, k:rb)) # Slice corollary upper side
                    end

                    proms[i] = abs(x[m[i]] - ($extrem)(coalesce(lref, rref), coalesce(rref, lref)))
                end
            end

            if isnothing(minprom)
                return (m, proms)
            else
                matched = findall(x -> x >= minprom, proms)
                return (m[matched], proms[matched])
            end
        end

    end
end

@doc """
    peakprom([Maxima()], x[, w=1]; minprom, strictbounds) => (idxs, proms)

Find the indices of local maxima and their prominences in `x`. `w` sets the minimum allowed
distance between maxima. `minprom` sets the minimum prominence of returned maxima.

Peak prominence is calculated as the distance between the current maxima and the larger of
the minimums of the left and right bounding intervals. Bounding intervals extend from the
next/previous index from the current maxima to the first element larger than or equal to
the current maxima, or the end of the signal, whichever comes first.

When bounding intervals contain `NaN` or `missing`, the reported prominence for that peak
will be `NaN` or `missing`, respectively, if `strictbounds = true`; if `strictbounds = false`,
the reference level for the contaminated bounding interval will be the minimum non-`NaN` or
`missing` value.

# Examples
```jldoctest
julia> x = rand(1000);

julia> ma, pa = peakprom(x);

julia> mi, pi = peakprom(Minima(), -x);

julia> @assert (mi == ma) && (pa == pi)

```

See also: [`argmaxima`](@ref), [`findmaxima`](@ref)
"""
peakprom(::Maxima, x, w)

peakprom(x, w=1; kwargs...) = peakprom(Maxima(), x, w; kwargs...)

@doc """
    peakprom(Minima(), x[, w=1]; minprom, strictbounds)

Find the indices of local minima and their prominences in `x`. `w` sets the minimum allowed
distance between minima. `minprom` sets the minimum prominence of returned minima.

Peak prominence is calculated as the distance between the current minima and the larger of
the minimums of the left and right bounding intervals. Bounding intervals extend from the
next/previous index from the current minima to the first element larger than or equal to
the current minima, or the end of the signal, whichever comes first.

When bounding intervals contain `NaN` or `missing`, the reported prominence for that peak
will be `NaN` or `missing`, respectively, if `strictbounds = true`; if `strictbounds = false`,
the reference level for the contaminated bounding interval will be the maximum non-`NaN` or
`missing` value.

See also: [`argminima`](@ref), [`findminima`](@ref)
"""
peakprom(::Minima, x, w)

