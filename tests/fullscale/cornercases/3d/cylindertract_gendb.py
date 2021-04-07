#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/cornercases/3d/cylindertract_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the pressurized cylinder test.

import math
import numpy


class GenerateDB(object):
    """
    Python script to generate spatial database with displacement
    boundary conditions for the pressurized cylinder test.
    """

    def run(self):
        """
        Generate the database.
        """
        # Domain
        a = 2000.0
        b = 5000.0
        r = numpy.linspace(a, b, num=5, dtype=numpy.float64)
        theta = numpy.linspace(0.0, 0.5*math.pi, num = 10, dtype=numpy.float64)
        z = numpy.linspace(-1.0e+3, 1.0e+3, num=3, dtype=numpy.float64)
        r3, t3, z3 = numpy.meshgrid(r, theta, z)
        ct = numpy.cos(t3)
        st = numpy.sin(t3)
        x3 = r3*ct
        y3 = r3*st
        npts = x3.ravel().shape[0]
        
        xyz = numpy.zeros((npts, 3), dtype=numpy.float64)
        xyz[:, 0] = x3.ravel()
        xyz[:, 1] = y3.ravel()
        xyz[:, 2] = z3.ravel()

        from cylindertract_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xyz)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
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

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("cylindertract_bc.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
