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

## @file unittests/pytests/meshio/TestOutputSolnSubset.py

## @brief Unit testing of Python OutputSolnSubset object.

import unittest

from pylith.meshio.OutputSolnSubset import OutputSolnSubset

# ----------------------------------------------------------------------
class TestOutputSolnSubset(unittest.TestCase):
  """
  Unit testing of Python OutputSolnSubset object.
  """

  def setUp(self):
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    iohandler = MeshIOAscii()
    filename = "data/twohex8.txt"
    
    from spatialdata.geocoords.CSCart import CSCart
    iohandler.filename = filename
    iohandler.coordsys = CSCart()
    mesh = iohandler.read(debug=False, interpolate=False)

    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)

    name = "solution"
    fields.addReal(name)
    fields.setFiberDimension(name, mesh.dimension())
    fields.allocate(name)
    fields.solutionField(name)

    self.mesh = mesh
    self.fields = fields
    return


  def test_constructor(self):
    """
    Test constructor.
    """
    output = OutputSolnSubset()
    output._configure()
    output.writer._configure()
    return


  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "label"
    output.preinitialize()
    
    self.assertEqual(output, output.dataProvider)
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "label"
    output.preinitialize()

    output.vertexDataFields = ["displacements"]
    output.verifyConfiguration()
    return
  
  
  def test_initialize(self):
    """
    Test initialize().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "2"
    output.writer._configure()
    output.writer.filename = "test.vtk"

    output.preinitialize()
    output.initialize(self.mesh)
    self.assertNotEqual(None, output.cppHandle)
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "2"
    output.writer._configure()
    output.writer.filename = "test.vtk"

    output.preinitialize()
    output.initialize(self.mesh)

    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.close()
    return


  def test_writeInfo(self):
    """
    Test writeInfo().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "2"
    output.writer._configure()
    output.writer.filename = "output_sub.vtk"

    output.preinitialize()
    output.initialize(self.mesh)
    
    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.writeInfo()
    output.close()
    return


  def test_writeData(self):
    """
    Test writeData().
    """
    output = OutputSolnSubset()
    output._configure()
    output.label = "2"
    output.writer._configure()
    output.writer.filename = "outputsub.vtk"
    output.writer.timeFormat = "%3.1f"
    output.vertexDataFields = ["displacements"]

    output.preinitialize()
    output.initialize(self.mesh)

    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.writeData(2.0*s, self.fields)
    output.close()
    return


# End of file 
