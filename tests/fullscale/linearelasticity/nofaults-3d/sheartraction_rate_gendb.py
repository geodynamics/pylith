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
# @file tests/fullscale/linearelasticity/nofaults-3d/sheartraction_rate_gendb.py
#
# @brief Python script to generate spatial database with Dirichlet
# boundary conditions for the time-dependent shear test. The traction
# boundary conditions use UniformDB in the .cfg file.

import numpy

from pythia.pyre.units.time import year


class GenerateDB(object):
    """Python object to generate spatial database with Dirichlet
    boundary conditions for the tim-dependent shear test.
    """

    def __init__(self):
        """Constructor.
        """
        return

    def run(self):
        """Generate the database.
        """
        # Domain
        x = numpy.arange(-1.0e+4, 1.01e+4, 5.0e+3)
        y = numpy.arange(-1.0e+4, 1.01e+4, 5.0e+3)
        z = numpy.array([0])
        x3, y3, z3 = numpy.meshgrid(x, y, z)
        nptsX = x.shape[0]
        nptsY = y.shape[0]
        nptsZ = z.shape[0]

        xyz = numpy.zeros((nptsX * nptsY * nptsZ, 3), dtype=numpy.float64)
        xyz[:, 0] = x3.ravel()
        xyz[:, 1] = y3.ravel()
        xyz[:, 2] = z3.ravel()

        from sheartraction_rate_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.bc_initial_displacement(xyz)
        velocity_time = soln.bc_rate_time(xyz) / year.value
        velocity = soln.bc_velocity(xyz) * year.value

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            "x": x,
            "y": y,
            "z": z,
            'points': xyz,
            'coordsys': cs,
            'data_dim': 2,
            'values': [
                {'name': "initial_amplitude_x",
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
                 'data': velocity_time[0, :, 0].ravel()},
            ]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "sheartraction_rate_disp.spatialdb"
        io._configure()
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
