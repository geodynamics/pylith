#!/usr/bin/env python
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
# @file tests/fullscale/linearelasticity/nofaults-3d/gravity_refstate_gendb.py
#
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
        z = numpy.arange(-9000.0, 0.1, 9000.0)
        x = numpy.zeros(z.shape)
        y = numpy.zeros(z.shape)
        npts = z.shape[0]

        xyz = numpy.vstack((x, y, z)).transpose()

        from gravity_refstate_soln import AnalyticalSoln
        from gravity_refstate_soln import p_density, p_vs, p_vp
        soln = AnalyticalSoln()
        stress = soln.stress(xyz)
        strain = soln.strain(xyz)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            'points': xyz,
            'coordsys': cs,
            'data_dim': 1,
            'values': [
                {
                    'name': "density",
                    'units': "kg/m**3",
                    'data': p_density * numpy.ones((npts,)),
                }, {
                    'name': "vs",
                    'units': "m/s",
                    'data': p_vs * numpy.ones((npts,)),
                }, {
                    'name': "vp",
                    'units': "m/s",
                    'data': p_vp * numpy.ones((npts,)),
                }, {
                    'name': "reference_stress_xx",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 0]),
                }, {
                    'name': "reference_stress_yy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 1]),
                }, {
                    'name': "reference_stress_zz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 2]),
                }, {
                    'name': "reference_stress_yz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 4]),
                }, {
                    'name': "reference_stress_xz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 5]),
                }, {
                    'name': "reference_stress_xy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0,:, 3]),
                }, {
                    'name': "reference_strain_xx",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 0]),
                }, {
                    'name': "reference_strain_yy",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 1]),
                }, {
                    'name': "reference_strain_zz",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 2]),
                }, {
                    'name': "reference_strain_xy",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 3]),
                }, {
                    'name': "reference_strain_yz",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 4]),
                }, {
                    'name': "reference_strain_xz",
                    'units': "none",
                    'data': numpy.ravel(strain[0,:, 5]),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("gravity_refstate_matfields.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
