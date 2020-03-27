export maxima,
       minima

for (funcname, comp) in ((:maxima, :<),
                         (:minima, :>))

    @eval begin
        function ($funcname)(x::AbstractVector,
                        w::Integer=1,
                        includebounds::Bool=false)
            w > 0 || throw(ArgumentError("window cannot be negative"))
            xlen = length(x)
            idxs = Int[]

            # There can't be more than one peak every `w` elements
            maxN = includebounds ? cld(xlen,2*w) : fld(xlen,2*w)
            N = 0

            ww = [ collect(w:-1:1); collect(-w:-1) ]

            sizehint!(idxs,maxN)

            i = findfirst(!ismissing, x)

            if includebounds
                @inbounds while i <= w
                    peak = true
                    for j in max(i-w, 1):(i+w)
                        i === j && continue
                        if coalesce(($comp)(x[i], x[j]), true) || x[i] === x[j]
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    if ismissing(x[i])
                        i = something(findnext(!ismissing, x, i+1), xlen+1)
                    else
                        i += 1
                    end
                end
            end

            i = something(findnext(!ismissing, x, 1+w), xlen+1)
            maxI = xlen-w

            @inbounds while i <= maxI
                peak = true
                for j in ww # For all elements within the window
                    if coalesce(($comp)(x[i], x[i-j]), true) || x[i] === x[i-j]
                        peak &= false
                        break # No reason to continue checking if the rest of the elements qualify
                    end
                end

                if peak
                    push!(idxs,i)
                    N += 1
                    i += w # There can't be another peak for at least `w` more elements
                end
                if ismissing(x[i])
                    i = something(findnext(!ismissing, x, i+1), xlen+1)
                else
                    i += 1
                end
            end

            if includebounds
                @inbounds while i <= xlen
                    peak = true
                    for j in (i-w):min((i+w),xlen)
                        i === j && continue
                        if coalesce(($comp)(x[i], x[j]), true) || x[i] === x[j]
                            peak &= false
                            break # No reason to continue checking if the rest of the elements qualify
                        end
                    end

                    if peak
                        push!(idxs,i)
                        N += 1
                        i += w # There can't be another peak for at least `w` more elements
                    end
                    if i < xlen && ismissing(x[i])
                        i = something(findnext(!ismissing, x, i+1), xlen+1)
                    else
                        i += 1
                    end
                end
            end

            resize!(idxs,N)
            sizehint!(idxs,N)

            return idxs
        end
    end
end

@doc """
    maxima(x[, w=1, includebounds=false])

Find the local maxima of `x`.

`w` sets the minimum allowed distance between maxima. If `includebounds` is `true`, maxima
are allowed to be less than `w` away from the bounds of `x`, otherwise, maxima may not be
located any closer to the ends of `x` than `w`.
"""
maxima

@doc """
    minima(x[, w=1, includebounds=false])

Find the local minima of `x`.

`w` sets the minimum allowed distance between minima. If `includebounds` is `true`, minima
are allowed to be less than `w` away from the bounds of `x`, otherwise, minima may not be
located any closer to the ends of `x` than `w`.
"""
minima

