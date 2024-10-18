#!/usr/bin/env nemesis

from matplotlib import pyplot
import numpy

from pythia.pyre.units.length import km

FILENAME = "output/slip_variable.txt"

data = numpy.loadtxt(FILENAME)
x = data[:,0] / km.value
slip = data[:,1]

figure = pyplot.figure(figsize=(6.5, 4.0))
ax = figure.add_axes(rect=(0.09, 0.12, 0.88, 0.86))
ax.spines[["top", "right"]].set_visible(False)
ax.plot(x, slip)
ax.set_xlabel("Distance along strike, km")
ax.set_ylabel("Left-lateral slip, m")

ax.set_xlim(numpy.min(x), numpy.max(x))
ax.set_ylim(0, 1.1*numpy.max(slip))

pyplot.savefig("step04-slip.pdf")
print(f"Mean slip: {numpy.mean(data[:,1]):.2f}")
