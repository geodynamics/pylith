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
# @file tests/fullscale/cornercases/3d/axialdisp_rotated_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the axial displacement test in rotated coordinates.

import math
import numpy


class GenerateDB(object):
    """
    Python object to generate spatial database with displacement
    boundary conditions for the axial displacement test in rotated coordinates.
    """

    def run(self):
        """
        Generate the database.
        """
        # Domain
        x = numpy.arange(-4.0e+3, 4.01e+3, 1.0e+3)
        y = numpy.arange(-4.0e+3, 4.01e+3, 1.0e+3)
        z = numpy.arange(-8.0e+3, 0.01e+3, 1.0e+3)
        x3, y3, z3 = numpy.meshgrid(x, y, z)
        nptsX = x.shape[0]
        nptsY = y.shape[0]
        nptsZ = z.shape[0]

        xyz = numpy.zeros((nptsX * nptsY * nptsZ, 3), dtype=numpy.float64)
        xyz[:, 0] = x3.ravel()
        xyz[:, 1] = y3.ravel()
        xyz[:, 2] = z3.ravel()

        # Rotation angles.
        alpha = numpy.radians(30.0) # Rotation about x axis.
        beta = numpy.radians(30.0) # Rotation about y axis.

        # Rotation matrix.
        ca = math.cos(alpha)
        sa = math.sin(alpha)
        cb = math.cos(beta)
        sb = math.sin(beta)
        rotMatA = numpy.array([[ ca, -sa, 0.0],
                               [ sa,  ca, 0.0],
                               [0.0, 0.0, 1.0]], dtype=numpy.float64)
        rotMatB = numpy.array([[ cb, 0.0,  sb],
                               [0.0, 1.0, 0.0],
                               [-sb, 0.0,  cb]], dtype=numpy.float64)
        rotMat = numpy.dot(rotMatA, rotMatB)

        # Rotated coordinates.
        xyzRotated = numpy.dot(xyz, rotMat.transpose())

        from axialdisp_rotated_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xyzRotated)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            "points": xyzRotated,
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
        io = createWriter("axialdisp_rotated_bc.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
