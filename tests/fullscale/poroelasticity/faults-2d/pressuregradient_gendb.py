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
# @file tests/fullscale/linearelasticity/nofaults-2d/axialdisp_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the axial displacement test.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with displacement
    boundary conditions for the axial displacement test.
    """

    def run(self):
        """Generate the databases."""
        # Domain
        x = numpy.arange(0.0, 8000.1, 1000.0)
        y = 0.0
        npts = x.shape[0]

        xy = numpy.zeros((npts, 2), dtype=numpy.float64)
        xy[:, 0] = x
        xy[:, 1] = y

        from pressuregradient_soln import AnalyticalSoln

        soln = AnalyticalSoln()
        disp = soln.displacement(xy)
        press = soln.pressure(xy)

        from spatialdata.geocoords.CSCart import CSCart

        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()

        disp = {
            "points": xy,
            "coordsys": cs,
            "data_dim": 1,
            "values": [
                {
                    "name": "initial_amplitude_x",
                    "units": "m",
                    "data": numpy.ravel(disp[0, :, 0]),
                },
                {
                    "name": "initial_amplitude_y",
                    "units": "m",
                    "data": numpy.ravel(disp[0, :, 1]),
                },
            ],
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter

        io = createWriter("pressuregradient_disp.spatialdb")
        io.write(disp)

        pressure = {
            "points": xy,
            "coordsys": cs,
            "data_dim": 1,
            "values": [
                {
                    "name": "initial_amplitude",
                    "units": "Pa",
                    "data": numpy.ravel(press[0, :, 0]),
                }
            ],
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter

        io = createWriter("pressuregradient_press.spatialdb")
        io.write(pressure)


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
