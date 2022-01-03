#!/usr/bin/env nemesis

import h5py
import matplotlib.pyplot as pyplot

fullpath = "/Users/baagaard/scratch/build/clang-13.0/cig/pylith-debug/tests/manual/petsc/sliderblock/sliderblock.h5"

h5 = h5py.File(fullpath, "r")
t = h5["time"][:].squeeze()
soln = h5["solution"][:]
u = soln[:, 0:4, :]
v = soln[:, 4:8, :]
l = soln[:, -1, :]

x = [0, 1, 2, 3]

for i in range(4):
    pyplot.plot(x[i] + u[:, i, :].squeeze(), t)
pyplot.show()

pyplot.plot(t, (u[:, 2, :]-u[:, 1, :]).squeeze())
pyplot.show()

pyplot.plot(t, (v[:, 2, :]-v[:, 1, :]).squeeze())
pyplot.show()

pyplot.plot(t, l.squeeze())
pyplot.show()
