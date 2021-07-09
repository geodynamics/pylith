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
# @file tests/fullscale/viscoelasticity/nofaults-3d/axialstrainrate_genmaxwell_gendb.py
#
# @brief Python script to generate spatial database with velocity
# boundary conditions for the axial strain test.

import numpy
from pythia.pyre.units.time import year


class GenerateDB(object):
    """Python object to generate spatial database with velocity
    boundary conditions for the axial strain test.
    """

    def __init__(self):
        """Constructor.
        """
        return

    def run(self):
        """Generate the database.
        """
        # Domain
        x = numpy.arange(-6000.0, 6000.1, 2000.0)
        y = numpy.arange(-6000.0, 6000.1, 2000.0)
        z = numpy.arange(-10000.0, 2000.1, 2000.0)
        (x3, y3, z3) = numpy.meshgrid(x, y, z)
        nptsX = x.shape[0]
        nptsY = y.shape[0]
        nptsZ = z.shape[0]

        xyz = numpy.zeros((nptsX * nptsY * nptsZ, 3), dtype=numpy.float64)
        xyz[:, 0] = x3.ravel()
        xyz[:, 1] = y3.ravel()
        xyz[:, 2] = z3.ravel()

        from axialstrainrate_genmaxwell_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.initial_displacement(xyz)
        velocity_time = soln.bc_rate_time(xyz)/year.value
        velocity = soln.bc_velocity(xyz)*year.value

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            'x': x,
            'y': y,
            'z': z,
            'points': xyz,
            'coordsys': cs,
            'data_dim': 3,
            'values': [{'name': "initial_amplitude_x",
                        'units': "m",
                        'data': disp[0, :, 0].ravel()},
                       {'name': "initial_amplitude_y",
                        'units': "m",
                        'data': disp[0, :, 1].ravel()},
                       {'name': "initial_amplitude_z",
                        'units': "m",
                        'data': disp[0, :, 2].ravel()},
                       {'name': "rate_amplitude_x",
                        'units': "m/year",
                        'data': velocity[0, :, 0].ravel()},
                       {'name': "rate_amplitude_y",
                        'units': "m/year",
                        'data': velocity[0, :, 1].ravel()},
                       {'name': "rate_amplitude_z",
                        'units': "m/year",
                        'data': velocity[0, :, 2].ravel()},
                       {'name': "rate_start_time",
                        'units': "year",
                        'data': velocity_time[0, :, 0].ravel()}
                       ]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "axialstrainrate_genmaxwell_bc.spatialdb"
        io._configure()
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
