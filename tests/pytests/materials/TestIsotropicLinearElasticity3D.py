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

# @file tests/pytests/materials/TestIsotropicLinearElasticity3D.py

# @brief Unit testing of TestIsotropicLinearElasticity3D object.

import unittest

from pylith.materials.IsotropicLinearElasticity3D import IsotropicLinearElasticity3D
from pylith.tests.UnitTestApp import configureSubcomponents


class TestIsotropicLinearElasticity3D(unittest.TestCase):
    """
    Unit testing of TestIsotropicLinearElasticity3D object.
    """

    def setUp(self):
        """
        Setup test subject.
        """
        self.material = IsotropicLinearElasticity3D()
        configureSubcomponents(self.material)
        return

    def test_constructor(self):
        """
        Test constructor.
        """
        self.assertEqual(0, self.material.materialId)
        return

    def tes_preinitialize(self):
        """
        Test preinitialize(). Set inventory and verify values from C++ object.
        """
        from pylith.topology.Mesh import Mesh
        mesh = Mesh()

        materialId = 4
        label = "material ABC"
        useInertia = True
        useBodyForce = False
        useReferenceState = True

        self.material.inventory.materialId = materialId
        self.material.inventory.label = label
        self.material.inventory.useInertia = useInertia
        self.material.inventory.useBodyForce = useBodyForce
        self.material.inventory.useReferenceState = useReferenceState
        self.material.preinitialize(mesh)

        from pylith.materials.IsotropicLinearElasticity3D import ModuleMaterial
        self.assertEqual(materialId, ModuleMaterial.id(self.material))
        self.assertEqual(label, ModuleMaterial.label(self.material))
        self.assertEqual(useInertia, ModuleMaterial.useInertia(self.material))
        self.assertEqual(useBodyForce, ModuleMaterial.useBodyForce(self.material))
        self.assertEqual(useReferenceState, ModuleMaterial.useReferenceState(self.material))
        return

    def test_factory(self):
        """
        Test factory method.
        """
        from pylith.materials.IsotropicLinearElasticity3D import material
        self.assertTrue(isinstance(material(), IsotropicLinearElasticity3D))
        return


# End of file
