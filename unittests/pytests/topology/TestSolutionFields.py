#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
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
    self.mesh = importer.read(normalizer, debug=False, interpolate=False)
    
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
    fields.add("field A");
    fields.add("field B");
    fields.add("field C");
    
    fields.solutionName("field B")
    return


  def test_solution(self):
    """
    Test solution().
    """
    fields = self.fields
    fields.add("field A");
    fields.add("field B");
    fields.add("field C");
    
    fields.solutionName("field B")
    solution = self.fields.solution()
    return


  def test_createHistory(self):
    """
    Test createHistory().
    """
    fields = self.fields
    fields.add("field A");
    fields.add("field B");
    fields.add("field C");
    
    fields.createHistory(["field B", "field A", "field C"])
    return


  def test_shiftHistory(self):
    """
    Test shiftHistory().
    """
    fields = self.fields
    fields.add("field A");
    fields.add("field B");
    fields.add("field C");
    
    fields.createHistory(["field B", "field A", "field C"])
    fields.shiftHistory()
    return


# End of file 
