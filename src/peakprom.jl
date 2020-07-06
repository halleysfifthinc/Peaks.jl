using UnsafeArrays

export Maxima,
       Minima

export peakprom

struct Maxima; end
struct Minima; end

for (extrema, comps, Extrema) in ((maxima, (:>, minimum, max), Maxima),
                       (minima, (:<, maximum, min), Minima))

    comp, argcomp1, argcomp2 = comps

    @eval begin
        function peakprom(x::AbstractVector{T}, ::$Extrema, w::Integer=1, minprom::T=zero(T), strictbounds::Bool=true) where T
            m = ($extrema)(x, w, strictbounds)
            # the peaks should always be compared to the peaks within the boundary
            mcomp = strictbounds ? ($extrema)(x, w, false) : m

            M = lastindex(m)
            Mcomp = lastindex(mcomp)
            lbegin = firstindex(x)
            lend = lastindex(x)
            proms = similar(x,M)

            @inbounds for i in eachindex(m)
                lb, rb = lbegin, lend
                lcan, rcan = zero(T), zero(T)

                ci = findnext(x -> x == m[i], mcomp, i)

                # Find left bound
                for ii in (ci - 1):-1:1
                    if ($comp)(x[mcomp[ii]], x[m[i]])
                        lb = mcomp[ii]
                        break
                    end
                end

                # Find right bound
                for ii in (ci + 1):Mcomp
                    if ($comp)(x[mcomp[ii]], x[m[i]])
                        rb = mcomp[ii]
                        break
                    end
                end

                if strictbounds
                    lcan = m[i] == lbegin ? x[m[i]] : ($argcomp1)(uview(x, lb:(m[i] - 1))) # Slice from lower bound to the index prior to the current extrema
                    rcan = m[i] == lend ? x[m[i]] : ($argcomp1)(uview(x, (m[i] + 1):rb)) # Slice corollary upper side
                else
                    lcan = m[i] == lbegin ? x[m[i]] : ($argcomp1)(filter(!isnan, (skipmissing(uview(x, lb:(m[i] - 1)))))) # Slice from lower bound to the index prior to the current extrema, ignoring NaNs and missings
                    rcan = m[i] == lend ? x[m[i]] : ($argcomp1)(filter(!isnan, (skipmissing(uview(x, (m[i] + 1):rb))))) # Slice corollary upper side
                end
                proms[i] = abs(x[m[i]] - ($argcomp2)(lcan, rcan))
            end

            if minprom != zero(T)
                matched = findall(x -> x >= minprom, proms)
                return (m[matched], proms[matched])
            else
                return (m, proms)
            end
        end

    end
end

@doc """
    peakprom(x, ::Maxima[, w=1, minprom=0])

Find all local maxima and maxima prominences in `x` matching the conditions `w` and `minprom`.
`w` sets the minimum allowed distance between maxima. `minprom` sets the minimum prominence
(inclusive) of returned maxima.

Peak prominence is calculated as the distance between the current maxima and the highest of
the minimums of the lower and upper bounds. Bounds extend from the next index from the
current maxima to the next maxima higher than the current maxima, or the end of the signal,
which ever comes first.

# Examples
```jldoctest
julia> x = rand(1000);

julia> ma, pa = peakprom(x, Maxima());

julia> mi, pi = peakprom(-x, Minima());

julia> @assert (mi == ma) && (pa == pi)

```

See also: [`maxima`](@ref), [`minima`](@ref)
"""
peakprom(x, ::Maxima)

@doc """
    peakprom(x, ::Minima[, w=1, minprom=0])

Find all local minima and minima prominences in `x` matching the conditions `w` and `minprom`.
`w` sets the minimum allowed distance between minima. `minprom` sets the minimum prominence
(inclusive) of returned minima.

Peak prominence is calculated as the distance between the current minima and the lowest of
the maximums of the lower and upper bounds. Bounds extend from the next index from the
current minima to the next minima lower than the current minima, or the end of the signal,
which ever comes first.

See also: [`maxima`](@ref), [`minima`](@ref)
"""
peakprom(x, ::Minima)

