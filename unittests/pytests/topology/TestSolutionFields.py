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

## @file unittests/pytests/topology/TestSolutionFields.py

## @brief Unit testing of SolutionFields object.

import unittest

from pylith.topology.SolutionFields import SolutionFields

# ----------------------------------------------------------------------
class TestSolutionFields(unittest.TestCase):
  """
  Unit testing of SolutionFields object.
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
    
    self.fields = SolutionFields(self.mesh)
    return


  def test_constructor(self):
    """
    Test constructor.
    """
    return


  def test_solutionName(self):
    """
    Test mesh().
    """
    fields = self.fields
    fields.add("field A", "A");
    fields.add("field B", "B");
    fields.add("field C", "C");
    
    fields.solutionName("field B")
    return


  def test_solution(self):
    """
    Test solution().
    """
    fields = self.fields
    fields.add("field A", "A");
    fields.add("field B", "B");
    fields.add("field C", "C");
    
    fields.solutionName("field B")
    solution = self.fields.solution()
    return


  def test_fieldAdd(self):
    """
    Test fieldAdd().
    """
    fields = self.fields
    fields.add("field A", "A");
    fields.add("field B", "B");

    helper_fieldAdd(fields)
    fieldA = fields.get("field A")
    fieldB = fields.get("field B")
    fieldA.allocate()
    fieldB.allocate()
    fieldA.copy(fieldB)
    return


def helper_fieldAdd(fields):
  fieldA = fields.get("field A")
  fieldB = fields.get("field B")
  fieldA.allocate()
  fieldB.allocate()

  fieldA.add(fieldB)
  return


# End of file 
