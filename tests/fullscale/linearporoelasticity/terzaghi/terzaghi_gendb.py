#!/usr/bin/env python
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
# @file tests/fullscale/linearporoelasticity/terzaghi/terzaghi_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the terzaghi test.

import numpy



class GenerateDB(object):
    """
    Python object to generate spatial database with initial conditions
    for the terzaghi poroelastic test.
    """

    def run(self):
        """
        Generate the database.
        """
        # Domain
        x = numpy.arange(0.0, 1.1, 0.1)
        y = numpy.arange(0.0, 10.1, 1.0)
        npts = x.shape[0]

        xx = x * numpy.ones((npts, 1), dtype=numpy.float64)
        yy = y * numpy.ones((npts, 1), dtype=numpy.float64)
        xy = numpy.zeros((npts**2, 2), dtype=numpy.float64)
        xy[:, 0] = numpy.ravel(xx)
        xy[:, 1] = numpy.ravel(numpy.transpose(yy))

        from terzaghi_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.initial_displacement(xy)
        pres = soln.initial_pressure(xy)
        trace_strain = soln.initial_trace_strain(xy)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {
            'x' : x,
            'y' : y,
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [{'name': "displacement_x",
                        'units': "m",
                        'data': numpy.ravel(disp[0,:, 0])},
                       {'name': "displacement_y",
                        'units': "m",
                        'data': numpy.ravel(disp[0,:, 1])},
                       {'name': "pressure",
                        'units': "Pa",
                        'data': numpy.ravel(pres[0,:])},
                       {'name': "trace_strain",
                        'units': "none",
                        'data': numpy.ravel(trace_strain[0,:])}]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "terzaghi_ic.spatialdb"
        io._configure()
        io.write(data)

        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
