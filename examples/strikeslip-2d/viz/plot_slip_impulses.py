#!/usr/bin/env nemesis

import pathlib

from matplotlib import pyplot

import numpy
import h5py


OUTPUT_DIR = pathlib.Path("output")
SIM_NAME = "step05_greensfns"

filename = f"{SIM_NAME}-fault.h5"

h5 = h5py.File(OUTPUT_DIR / filename)
y = h5["geometry/vertices"][:,1]
slip = h5["vertex_fields/slip"][:,:,1]
n_impulses = slip.shape[0]

i_sort = numpy.argsort(y)
y = y[i_sort] / 1.0e+3
slip = slip[:,i_sort]

figure = pyplot.figure(figsize=(6.5, 4.0))
ax = figure.add_axes(rect=(0.09, 0.13, 0.9, 0.83))
ax.spines[["top", "right"]].set_visible(False)
for i_impulse in range(n_impulses):
    ax.plot(y, slip[i_impulse,:])
ax.set_xlabel("Distance along strike, km")
ax.set_ylabel("Impulse amplitude, m")
ax.autoscale(axis="both", enable=True, tight="True")

pyplot.savefig(f"{SIM_NAME}-impulses.pdf")
pyplot.close(figure)
