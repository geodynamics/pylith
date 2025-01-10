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
# @file tests/fullscale/poroelasticity/mandel_compaction/mandel_compaction_gendb.py
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

        from mandel_compaction_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.initial_displacement(xy)
        pres = soln.initial_pressure(xy)
        trace_strain = soln.initial_trace_strain(xy)
        vel = soln.zero_array_dim(xy)
        pres_t = soln.zero_array(xy)
        trace_strain_t = soln.zero_array(xy)

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
            'values': [{'name': "displacement_x",
                        'units': "m",
                        'data': numpy.ravel(disp[0, :, 0])},
                       {'name': "displacement_y",
                        'units': "m",
                        'data': numpy.ravel(disp[0, :, 1])},
                       {'name': "pressure",
                        'units': "Pa",
                        'data': numpy.ravel(pres[0, :])},
                       {'name': "trace_strain",
                        'units': "none",
                        'data': numpy.ravel(trace_strain[0, :])},
                       {'name': "velocity_x",
                        'units': "m/s",
                        'data': numpy.ravel(vel[0, :, 0])},
                       {'name': "velocity_y",
                        'units': "m/s",
                        'data': numpy.ravel(vel[0, :, 1])},
                       {'name': "pressure_t",
                        'units': "Pa/s",
                        'data': numpy.ravel(pres_t[0, :])},
                       {'name': "trace_strain_t",
                        'units': "1/s",
                        'data': numpy.ravel(trace_strain_t[0, :])}]}

        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        io = SimpleGridAscii()
        io.inventory.filename = "mandel_compaction_bc.spatialdb"
        io._configure()
        io.write(data)
        data["values"] = [
            {
                'name': "displacement_x",
                'units': "m",
                'data': numpy.ravel(disp[0, :, 0])
            },  {
                'name': "displacement_y",
                'units': "m",
                'data': numpy.ravel(disp[0, :, 1])
            },  {
                'name': "pressure",
                'units': "Pa",
                'data': numpy.ravel(pres[0, :])
            },  {
                'name': "trace_strain",
                'units': "none",
                'data': numpy.ravel(trace_strain[0, :])
            },  {
                'name': "velocity_x",
                'units': "m/s",
                'data': numpy.ravel(vel[0, :, 0])
            },  {
                'name': "velocity_y",
                'units': "m/s",
                'data': numpy.ravel(vel[0, :, 1])
            },  {
                'name': "pressure_t",
                'units': "Pa/s",
                'data': numpy.ravel(pres_t[0, :])
            },  {
                'name': "trace_strain_t",
                        'units': "1/s",
                        'data': numpy.ravel(trace_strain_t[0, :])
            }]
        io.inventory.filename = "mandel_compaction_ic.spatialdb"
        io._configure()
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
