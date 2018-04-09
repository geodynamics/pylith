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

# @file tests_auto/linearelasticity/nofaults/gravity_nodeform_gendb.py
##
# @brief Python script to generate spatial database with auxiliary
# fields for test with gravitational body forces and initial
# stress/strain but no displacement.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with auxiliary fields for
    test with gravitational body forces and initial stress/strain but no
    displacement.
    """

    def __init__(self):
        """
        Constructor.
        """
        return

    def run(self):
        """
        Generate the database.
        """
        # Domain
        y = numpy.arange(-4000.0, 4000.1, 8000.0)
        x = numpy.zeros(y.shape)

        xy = numpy.vstack((x,y)).transpose()

        from gravity_nodeform_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        stress = soln.stress(xy)
        strain = soln.strain(xy)
        matprops = soln.matprops(xy)
        

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 1,
            'values': [
                {
                    'name': "density",
                    'units': "kg/m**3",
                    'data': numpy.ravel(matprops["density"]),
                },{
                    'name': "vs",
                    'units': "m/s",
                    'data': numpy.ravel(matprops["vs"]),
                },{
                    'name': "vp",
                    'units': "m/s",
                    'data': numpy.ravel(matprops["vp"]),
                },{
                    'name': "reference_stress_xx",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 0]),
                },{
                    'name': "reference_stress_yy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 1]),
                },{
                    'name': "reference_stress_zz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 2]),
                },{
                    'name': "reference_stress_xy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 3]),
                },{
                    'name': "reference_strain_xx",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 0]),
                },{
                    'name': "reference_strain_yy",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 1]),
                },{
                    'name': "reference_strain_zz",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 2]),
                },{
                    'name': "reference_strain_xy",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 3]),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
        io = SimpleIOAscii()
        io.inventory.filename = "gravity_nodeform_matfields.spatialdb"
        io._configure()
        io.write(data)
        return

# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
