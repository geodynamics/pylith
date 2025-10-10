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
        import mandel_soln
        soln = mandel_soln.AnalyticalSoln()
        locs = numpy.array([[0, mandel_soln.y_max]])
        t_steps = numpy.arange(0.0, 100.0e+5, 1.0e+4)

        # Need to restore tsteps in mandel_soln for checking fields
        soln_tsteps = mandel_soln.tsteps
        mandel_soln.tsteps = t_steps
        displacement = soln.displacement(locs)

        from spatialdata.spatialdb.TimeHistoryIO import write
        write(t_steps, displacement[:,0,1], units="m", filename="mandel_disp.timedb")
        mandel_soln.tsteps = soln_tsteps


# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
