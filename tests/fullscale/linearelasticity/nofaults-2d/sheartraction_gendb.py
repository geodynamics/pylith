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
# @file tests/fullscale/linearelasticity/nofaults-2d/sheartraction_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the shear test. The traction boundary
# conditions use UniformDB in the .cfg file.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with displacement
    boundary conditions for the shear test.
    """

    def __init__(self):
        """Constructor.
        """
        return

    def run(self):
        """Generate the database.
        """
        # Domain
        x = numpy.arange(-4000.0, 4000.1, 1000.0)
        y = numpy.arange(-4000.0, 4000.1, 1000.0)
        npts = x.shape[0]

        xx = x * numpy.ones((npts, 1), dtype=numpy.float64)
        yy = y * numpy.ones((npts, 1), dtype=numpy.float64)
        xy = numpy.zeros((npts**2, 2), dtype=numpy.float64)
        xy[:, 0] = numpy.ravel(xx)
        xy[:, 1] = numpy.ravel(numpy.transpose(yy))

        from sheartraction_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xy)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {'points': xy,
                'coordsys': cs,
                'data_dim': 2,
                'values': [{'name': "initial_amplitude_x",
                            'units': "m",
                            'data': disp[0, :, 0].ravel()},
                           {'name': "initial_amplitude_y",
                            'units': "m",
                            'data': disp[0, :, 1].ravel()}]}

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("sheartraction_disp.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
