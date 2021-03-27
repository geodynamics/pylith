#!/usr/bin/env nemesis
"""
This script generates a plot showing slip or fault tractions.
"""

# The code requires the numpy, h5py, and matplotlib packages.
import matplotlib.pyplot as pyplot
import numpy
import h5py
import matplotlib
matplotlib.use('Agg')


# ----------------------------------------------------------------------
def calcDist(vertices):
    """Compute down-dip distance from the trench.
    """
    dist = numpy.zeros(vertices.shape[0])
    pt0 = vertices[:-1, :]
    pt1 = vertices[1:, :]
    dx = ((pt1[:, 0]-pt0[:, 0])**2 + (pt1[:, 1]-pt0[:, 1])**2)**0.5
    dist[1:] = numpy.cumsum(dx)
    return dist

# ----------------------------------------------------------------------


def getData(sim):
    """Read fault information from HDF5 file.
    """
    filename = "output/%s-fault-slabtop.h5" % sim
    h5 = h5py.File(filename, "r")
    vertices = h5['geometry/vertices'][:]
    slip = h5['vertex_fields/slip'][:, :, :]
    tstamps = h5["time"][:]
    h5.close()

    data = {
        "time": tstamps,
        "vertices": vertices,
        "slip": slip
    }
    return data

# ----------------------------------------------------------------------


def plot(sim):

    # Get fault data for simulation.
    data = getData(sim)

    # Create sort key corresponding to increasing depth.
    indices = numpy.argsort(data["vertices"][:, 1])[::-1]

    # Calculate down-dip distance from trench and get sorted data.
    #dist = calcDist(data["vertices"][indices,:])
    dist = -data["vertices"][indices, 1]
    slip = data["slip"][:, indices, :]

    figure = pyplot.figure(figsize=(5.0, 3.0), facecolor='white', dpi=150)
    figure.set_facecolor('white')

    axes = figure.add_axes([0.15, 0.15, 0.80, 0.82])

    for i, _ in enumerate(data["time"]):
        color = "blue"
        lw = 0.5
        if i % 10 == 0:
            color = "red"
            lw = 1.0
        axes.plot(-slip[i, :, 0], dist/1.0e+3, '-', color=color, linewidth=lw)
    axes.set_xlabel("Slip (m)")
    #axes.set_ylabel("Down-dip Dist. (km)")
    axes.set_ylabel("Depth (km)")
    axes.invert_yaxis()

    pyplot.show()
    pyplot.savefig("subduction2d_%s_slip.pdf" % sim)
    return


# ======================================================================
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim", default="step05")
    args = parser.parse_args()

    plot(args.sim)


# End of file
