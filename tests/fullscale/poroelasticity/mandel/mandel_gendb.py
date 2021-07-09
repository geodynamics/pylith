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
# @file tests/fullscale/poroelasticity/mandel/mandel_gendb.py
#
# @brief Python script to generate spatial database with displacement
# boundary conditions for the mandel test.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with displacement
    boundary conditions for the axial displacement test.
    """

    def run(self):
        """Generate the database.
        """
        # Domain
        x1 = numpy.arange(-0.1, 10.1, 0.1)
        y1 = numpy.arange(-0.1, 1.01, 0.1)
        x, y = numpy.meshgrid(x1, y1)

        xy = numpy.zeros((len(x1) * len(y1), 2), dtype=numpy.float64)
        xy[:, 0] = x.ravel()
        xy[:, 1] = y.ravel()

        from mandel_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.initial_displacement(xy)
        pres = soln.initial_pressure(xy)
        trace_strain = soln.initial_trace_strain(xy)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {
            'x': x1,
            'y': y1,
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [{'name': "initial_amplitude_x",
                        'units': "m",
                        'data': numpy.ravel(disp[0, :, 0])},
                       {'name': "initial_amplitude_y",
                        'units': "m",
                        'data': numpy.ravel(disp[0, :, 1])},
                       {'name': "initial_pressure",
                        'units': "Pa",
                        'data': numpy.ravel(pres[0, :])},
                       {'name': "initial_trace_strain",
                        'units': "none",
                        'data': numpy.ravel(trace_strain[0, :])}]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "mandel_bc.spatialdb"
        io._configure()
        io.write(data)
        data["values"] = [
            {
                'name': "displacement_x",
                'units': "m",
                'data': numpy.ravel(disp[0, :, 0])
            }, {
                'name': "displacement_y",
                'units': "m",
                'data': numpy.ravel(disp[0, :, 1])
            }, {
                'name': "pressure",
                'units': "Pa",
                'data': numpy.ravel(pres[0, :])
            }, {
                'name': "trace_strain",
                'units': "none",
                'data': numpy.ravel(trace_strain[0, :])
            }]
        io.inventory.filename = "mandel_ic.spatialdb"
        io._configure()
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
