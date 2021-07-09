#!/usr/bin/env nemesis
#
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
# @file tests/fullscale/linearelasticity/nofaults-3d/axialdisp_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the axial displacement test.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with displacement
    boundary conditions for the axial displacement test.
    """

    def run(self):
        """Generate the database.
        """
        # Domain
        x = numpy.arange(-1.0e+4, 1.01e+4, 5.0e+3)
        y = numpy.arange(-1.0e+4, 1.01e+4, 5.0e+3)
        z = numpy.arange(-1.0e+4, 0.01e+4, 5.0e+3)
        x3, y3, z3 = numpy.meshgrid(x, y, z)
        nptsX = x.shape[0]
        nptsY = y.shape[0]
        nptsZ = z.shape[0]

        xyz = numpy.zeros((nptsX * nptsY * nptsZ, 3), dtype=numpy.float64)
        xyz[:, 0] = x3.ravel()
        xyz[:, 1] = y3.ravel()
        xyz[:, 2] = z3.ravel()

        from axialdisp_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xyz)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            "x": x,
            "y": y,
            "z": z,
            "points": xyz,
            "coordsys": cs,
            "data_dim": 3,
            "values": [
                {"name": "initial_amplitude_x",
                 "units": "m",
                 "data": numpy.ravel(disp[0, :, 0])},
                {"name": "initial_amplitude_y",
                 "units": "m",
                 "data": numpy.ravel(disp[0, :, 1])},
                {"name": "initial_amplitude_z",
                 "units": "m",
                 "data": numpy.ravel(disp[0, :, 2])},
            ]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "axialdisp_bc.spatialdb"
        io._configure()
        io.write(data)

        data["values"][0]["name"] = "displacement_x"
        data["values"][1]["name"] = "displacement_y"
        data["values"][2]["name"] = "displacement_z"
        io.inventory.filename = "axialdisp_ic.spatialdb"
        io._configure()
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
