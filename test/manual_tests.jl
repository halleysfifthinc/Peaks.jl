# These tests are intended to be run manually and interactively 
# to investigate the current functioning of the package.

data = [1, 2, 3, 4, 5, 4, 3, 2, 1, 6, 1]
pks = findmaxima(data)
pks = peakproms!(pks)
pks = peakwidths!(pks)

## Below is code intended to generate docstring examples
data = [1, 5, 1, 3, 2]
pks = findmaxima(data)
pks = peakproms!(pks)
pks = peakwidths!(pks)

##
run_tests = true
if run_tests
    using Pkg
    pkg_path = joinpath(homedir(), ".julia", "dev", "Peaks")
    Pkg.activate("PeaksTestEnv"; shared=true)
    Pkg.develop(path=pkg_path)
    Pkg.add(["OffsetArrays", "Plots"])
    include(joinpath(pkg_path, "test", "runtests.jl"))
end

##
begin
    fs = 100
    T = 1/fs
    t = 0:T:6pi+T
    sint = sin.(t)

    begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, sint[sinpks]; strict=false, relheight=1)
        @assert isapprox.(widths, fill(pi*100, length(sinpks)), atol=.01)|>all

        sinpks = argminima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, abs.(sint[sinpks]); strict=false, relheight=1)
        @assert isapprox.(widths, fill(pi*100, length(sinpks)), atol=.01)|>all
    end

    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [1.])
    @assert widths == [1.]
    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [NaN])
    @assert widths[1] === NaN
    _, widths, _, _ = peakwidths([2], [0.,1.,0.], [missing])
    @assert widths[1] === missing
    _, widths, _, _ = peakwidths([2], [0.,1.,NaN], [1.]; strict=true)
    @assert widths[1] === NaN
    _, widths, _, _ = peakwidths([2], [0.,1.,0.,-1.], [1.]; strict=false)
    _, widthsnan, _, _ = peakwidths([2], [0.,1.,NaN,-1.], [1.]; strict=false)
    @assert widths == widthsnan

    begin
        sinpks = argmaxima(sint)
        _, widths, _, _ = peakwidths(sinpks, sint, ones(length(sinpks)); strict=true, relheight=1)
        @assert first(widths) === NaN
        _, widths, _, _ = peakwidths(sinpks, sint, ones(length(sinpks)); strict=false, relheight=1)
        @assert first(widths) !== NaN
    end

    begin
        sinpks = argmaxima(sint)
        _, proms = peakproms(sinpks, sint)

        @assert length(first(peakwidths(sinpks, sint, proms; minwidth=pi*75))) == 2
        @assert length(first(peakwidths(sinpks, sint, proms; maxwidth=pi*75))) == 1
    end

end
