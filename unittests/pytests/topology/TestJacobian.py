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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/topology/TestJacobian.py

## @brief Unit testing of Jacobian object.

import unittest

from pylith.topology.Jacobian import Jacobian


# ----------------------------------------------------------------------
class TestJacobian(unittest.TestCase):
  """
  Unit testing of Jacobian object.
  """

  def setUp(self):
    """
    Setup mesh and associated field.
    """
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()    

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    importer = MeshIOAscii()
    importer.inventory.filename = "data/tri3.mesh"
    importer.inventory.coordsys = cs
    importer._configure()
    self.mesh = importer.read(debug=False, interpolate=False)

    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(self.mesh)
    fields.add("disp t+dt", "displacement")
    fields.solutionName("disp t+dt")
    solution = fields.solution()
    solution.newSection(solution.VERTICES_FIELD, self.mesh.dimension())
    solution.allocate()
    solution.zero()

    self.fields = fields
    self.jacobian = Jacobian(solution)
    return


  def tearDown(self):
    self.jacobian.cleanup()
    return


  def test_constructor(self):
    """
    Test constructor.
    """
    # setUp() tests constructor with default type
    jacobianA = Jacobian(self.fields.solution(), "aij")
    jacobianB = Jacobian(self.fields.solution(), "baij")
    return


  def test_matrix(self):
    """
    Test matrix().

    :WARNING: This is not a complete test of matrix(). We do not
    verify the results.
    """
    matrix = self.jacobian.matrix()

    # No testing of result.
    return


  def test_assemble(self):
    """
    Test assemble().

    :WARNING: This is not a complete test of assemble(). We do not
    verify the results.
    """
    self.jacobian.assemble("flush_assembly")
    self.jacobian.assemble("final_assembly")

    # No testing of result.
    return


  def test_zero(self):
    """
    Test zero().

    :WARNING: This is not a complete test of zero(). We do not
    verify the results.
    """
    self.jacobian.zero()

    # No testing of result.
    return


  def test_view(self):
    """
    Test view().

    :WARNING: This is not a complete test of view(). We do not
    verify the results.
    """
    self.jacobian.assemble("final_assembly")
    self.jacobian.view()

    # No testing of result.
    return


  def test_write(self):
    """
    Test write().

    :WARNING: This is not a complete test of write(). We do not
    verify the results.
    """
    self.jacobian = Jacobian(self.fields.solution(), "aij")
    self.jacobian.assemble("final_assembly")

    self.jacobian.write("jacobian.mat", self.mesh.comm())

    # No testing of result.
    return


# End of file 
