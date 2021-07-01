#!/usr/bin/env nemesis
# -*- Python -*- (syntax highlighting)
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# This script creates a spatial database with initial stresses for
# Step08b and Step08c using the PyLith output from Step08a.


from spatialdata.spatialdb.SimpleIOAscii import createWriter
import numpy
import h5py
import sys

sys.path.append('../mesh')


def getCellCenters(vertices, cells):
    """Function to compute cell centers.
    """
    cellCoords = vertices[cells,:]
    cellCenters = numpy.mean(cellCoords, axis=1)
    return cellCenters


def generate(sim, fileRoot, materials):
    from coordsys import cs_mesh

    for material in materials:

        filenameH5 = "../output/%s-%s.h5" % (sim, material)
        filenameDB = "%s-%s.spatialdb" % (fileRoot, material)

        # Open HDF5 file and get coordinates, cells, and stress.
        h5 = h5py.File(filenameH5, "r")
        vertices = h5['geometry/vertices'][:]
        cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int)

        # Get stresses from final time step.
        stress = h5['cell_fields/stress'][-1,:,:]
        h5.close()

        # Compute coordinates of quadrature points.
        quadCoords = getCellCenters(vertices, cells)

        # Create writer for spatial database file
        values = [{'name': "stress-xx",
                   'units': "Pa",
                   'data': stress[:, 0]},
                  {'name': "stress-yy",
                   'units': "Pa",
                   'data': stress[:, 1]},
                  {'name': "stress-zz",
                   'units': "Pa",
                   'data': stress[:, 2]},
                  {'name': "stress-xy",
                   'units': "Pa",
                   'data': stress[:, 3]},
                  {'name': "stress-yz",
                   'units': "Pa",
                   'data': stress[:, 4]},
                  {'name': "stress-xz",
                   'units': "Pa",
                   'data': stress[:, 5]},
                  ]

        writer = createWriter(filenameDB)
        writer.write({
            'points': quadCoords,
            'coordsys': cs_mesh(),
            'data_dim': 3,
            'values': values
        })


# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim", default="step08a_gravity_refstate")
    parser.add_argument("--file-root", action="store", dest="fileRoot", default="mat_initial_stress_grav")
    parser.add_argument("--materials", action="store", dest="materials", default="crust,mantle,slab,wedge")
    args = parser.parse_args()

    materials = args.materials.split(",")
    generate(args.sim, args.fileRoot, materials)


# End of file
