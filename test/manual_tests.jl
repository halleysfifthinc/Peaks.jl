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