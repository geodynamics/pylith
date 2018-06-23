#!/usr/bin/env python
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

# @file unittests/pytests/materials/TestIsotropicLinearElasticityPlaneStrain.py

# @brief Unit testing of TestIsotropicLinearElasticityPlaneStrain object.

import unittest

from pylith.materials.IsotropicLinearElasticityPlaneStrain import IsotropicLinearElasticityPlaneStrain

# ----------------------------------------------------------------------


def configureSubcomponents(facility):
    """Configure subcomponents."""
    for component in facility.components():
        configureSubcomponents(component)
        component._configure()
    return


class TestIsotropicLinearElasticityPlaneStrain(unittest.TestCase):
    """
    Unit testing of TestIsotropicLinearElasticityPlaneStrain object.
    """

    def setUp(self):
        """
        Setup test subject.
        """
        self.material = IsotropicLinearElasticityPlaneStrain()
        configureSubcomponents(self.material)
        return

    def test_constructor(self):
        """
        Test constructor.
        """
        self.assertEqual(2, self.material.dimension())
        return

    def testPreinitialize(self):
        """
        Test preinitialize().
        """
        from pylith.topology.Mesh import Mesh
        mesh = Mesh()
        self.material.preinitialize(mesh)

        from pylith.materials.IsotropicLinearElasticityPlaneStrain import ModuleMaterial
        self.assertEqual(0, ModuleMaterial.id(self.material))
        self.assertEqual("", ModuleMaterial.label(self.material))
        self.assertEqual(False, ModuleMaterial.useInertia(self.material))
        return

    def test_factory(self):
        """
        Test factory method.
        """
        from pylith.materials.IsotropicLinearElasticityPlaneStrain import material
        m = material()
        return


# End of file
