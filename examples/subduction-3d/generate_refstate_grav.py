#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
# This script creates a spatial database with reference stresses for
# Step08c using the PyLith output from Step08b.


import numpy
import h5py

from spatialdata.spatialdb.SimpleIOAscii import createWriter

from coordsys import cs_mesh

# Elastic properties matching mat_*_elastic.spatialdb
elastic_properties = {
    "crust": {
        "density": 3000.0,
        "vs": 4.0,
        "vp": 7.0,
    },
    "wedge": {
        "density": 2400.0,
        "vs": 3.5,
        "vp": 6.2,
    },
    "slab": {
        "density": 3400.0,
        "vs": 4.5,
        "vp": 8.0,
    },
    "mantle": {
        "density": 3300.0,
        "vs": 4.2,
        "vp": 7.5,
    },
}


def getCellCenters(vertices, cells):
    """Function to compute cell centers."""
    cellCoords = vertices[cells, :]
    cellCenters = numpy.mean(cellCoords, axis=1)
    return cellCenters


def generate(sim, fileRoot, materials):

    for material in materials:

        filenameH5 = "output/%s-%s.h5" % (sim, material)
        filenameDB = "%s-%s.spatialdb" % (fileRoot, material)

        # Open HDF5 file and get coordinates of vertices and cells.
        h5 = h5py.File(filenameH5, "r")
        vertices = h5["geometry/vertices"][:]
        cells = numpy.array(h5["viz/topology/cells"][:], dtype=numpy.int64)

        # Get stresses from final time step.
        stress = h5["cell_fields/cauchy_stress"][-1, :, :]
        h5.close()

        # Compute coordinates of cells.
        cellCoords = getCellCenters(vertices, cells)

        # Create writer for spatial database file
        values = [
            {
                "name": "density",
                "units": "kg/m**3",
                "data": elastic_properties[material]["density"]
                * numpy.ones_like(stress[:, 0]),
            },
            {
                "name": "vs",
                "units": "km/s",
                "data": elastic_properties[material]["vs"]
                * numpy.ones_like(stress[:, 0]),
            },
            {
                "name": "vp",
                "units": "km/s",
                "data": elastic_properties[material]["vp"]
                * numpy.ones_like(stress[:, 0]),
            },
            {"name": "reference_stress_xx", "units": "Pa", "data": stress[:, 0]},
            {"name": "reference_stress_yy", "units": "Pa", "data": stress[:, 1]},
            {"name": "reference_stress_zz", "units": "Pa", "data": stress[:, 2]},
            {"name": "reference_stress_xy", "units": "Pa", "data": stress[:, 3]},
            {"name": "reference_stress_yz", "units": "Pa", "data": stress[:, 4]},
            {"name": "reference_stress_xz", "units": "Pa", "data": stress[:, 5]},
            {
                "name": "reference_strain_xx",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 0],
            },
            {
                "name": "reference_strain_yy",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 1],
            },
            {
                "name": "reference_strain_zz",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 2],
            },
            {
                "name": "reference_strain_xy",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 3],
            },
            {
                "name": "reference_strain_yz",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 4],
            },
            {
                "name": "reference_strain_xz",
                "units": "Pa",
                "data": numpy.zeros_like(stress)[:, 5],
            },
        ]

        writer = createWriter(filenameDB)
        writer.write(
            {
                "points": cellCoords,
                "coordsys": cs_mesh(),
                "data_dim": 3,
                "values": values,
            }
        )


# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sim", action="store", dest="sim", default="step08b_gravity_incompressible"
    )
    parser.add_argument(
        "--file-root",
        action="store",
        dest="fileRoot",
        default="mat_refstate_grav",
    )
    parser.add_argument(
        "--materials",
        action="store",
        dest="materials",
        default="crust,mantle,slab,wedge",
    )
    args = parser.parse_args()

    materials = args.materials.split(",")
    generate(args.sim, args.fileRoot, materials)


# End of file
