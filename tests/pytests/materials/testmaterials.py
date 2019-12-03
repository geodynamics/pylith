#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

# @file tests/materials/testmaterials.py

# @brief Python application for testing materials code.

from pylith.tests.UnitTestApp import UnitTestApp

import unittest


class TestApp(UnitTestApp):
    """
    Test application.
    """

    def __init__(self):
        """
        Constructor.
        """
        UnitTestApp.__init__(self)
        return

    def _suite(self):
        """
        Setup the test suite.
        """

        suite = unittest.TestSuite()

        # Plane strain
        from TestIsotropicLinearElasticityPlaneStrain import TestIsotropicLinearElasticityPlaneStrain
        suite.addTest(unittest.makeSuite(TestIsotropicLinearElasticityPlaneStrain))

        from TestIsotropicLinearMaxwellPlaneStrain import TestIsotropicLinearMaxwellPlaneStrain
        suite.addTest(unittest.makeSuite(TestIsotropicLinearMaxwellPlaneStrain))

        from TestIsotropicLinearGenMaxwellPlaneStrain import TestIsotropicLinearGenMaxwellPlaneStrain
        suite.addTest(unittest.makeSuite(TestIsotropicLinearGenMaxwellPlaneStrain))

        # 3D
        from TestIsotropicLinearElasticity3D import TestIsotropicLinearElasticity3D
        suite.addTest(unittest.makeSuite(TestIsotropicLinearElasticity3D))

        from TestIsotropicLinearMaxwell3D import TestIsotropicLinearMaxwell3D
        suite.addTest(unittest.makeSuite(TestIsotropicLinearMaxwell3D))

        from TestIsotropicLinearGenMaxwell3D import TestIsotropicLinearGenMaxwell3D
        suite.addTest(unittest.makeSuite(TestIsotropicLinearGenMaxwell3D))

        # Containers
        from TestHomogeneous import TestHomogeneous
        suite.addTest(unittest.makeSuite(TestHomogeneous))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestApp().run()


# End of file
