# These tests are intended to be run manually and interactively 
# to investigate the current functioning of the package.

data = [1, 2, 3, 4, 5, 4, 3, 2, 1, 6, 1]
pks = findmaxima(data)
pks = peakproms!(pks)
pks = peakwidths!(pks)