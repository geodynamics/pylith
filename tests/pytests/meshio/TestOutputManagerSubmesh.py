#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#

## @file tests/pytests/meshio/TestOutputManagerSubmesh.py

## @brief Unit testing of Python OutputManagerSubmesh object.

import unittest

from pylith.meshio.OutputManager import OutputManager

# ----------------------------------------------------------------------
class TestProvider(object):

  def __init__(self):
    """Constructor.
    """
    self.availableFields = \
        {'vertex': \
           {'info': ["vertex info"],
            'data': ["vertex data 1",
                     "vertex data 2"]},
         'cell': \
           {'info': ["cell info"],
            'data': ["cell data"]}}

    filename = "data/twohex8.txt"
    
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    iohandler = MeshIOAscii()
    iohandler.inventory.filename = filename
    from spatialdata.geocoords.CSCart import CSCart
    iohandler.inventory.coordsys = CSCart()
    iohandler._configure()

    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()
    mesh = iohandler.read(debug=False, interpolate=False)

    from pylith.topology.Mesh import Mesh
    submesh = Mesh(mesh=mesh, label="4")

    from pylith.topology.Fields import Fields
    fields = Fields(submesh)
    
    self.mesh = mesh
    self.submesh = submesh
    self.fields = fields
    return


  def getDataMesh(self):
    """Get mesh.
    """
    return (self.submesh, None, None)


  def getVertexField(self, name, fields=None):
    """Get vertex field.
    """
    from pylith.field.field.FieldBase import SCALAR, VECTOR, OTHER
    if name == "vertex info":
      fieldType = FieldBase.SCALAR
      fiberDim = 1
    elif name == "vertex data 1":
      fieldType = FieldBase.VECTOR
      fiberDim = self.mesh.dimension()
    elif name == "vertex data 2":
      fieldType = FieldBase.OTHER
      fiberDim = 5
    else:
      raise ValueError("Unknown field '%s'." % name)

    self.fields.add(name, name)
    field = self.fields.get(name)
    field.newSection(field.VERTICES_FIELD, fiberDim)
    field.allocate()
    field.vectorFieldType(fieldType)
    self.fields.setFiberDimension(name, fiberDim)
    self.fields.allocate(name)
    field = self.fields.getReal(name)
    return field


  def getCellField(self, name, fields=None):
    """Get cell field.
    """
    if name == "cell info":
      fieldType = FieldBase.SCALAR
      fiberDim = 1
    elif name == "cell data":
      fieldType = FieldBase.VECTOR
      fiberDim = self.mesh.dimension()
    else:
      raise ValueError("Unknown field '%s'." % name)

    self.fields.add(name, name)
    field = self.fields.get(name)
    field.newSection(field.CELLS_FIELD, fiberDim)
    field.allocate()
    field.vectorFieldType(fieldType)
    return field


# ----------------------------------------------------------------------
class TestOutputManagerSubmesh(unittest.TestCase):
  """Unit testing of Python MeshOutputManager object.
  """

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    self.normalizer = Nondimensional()
    self.normalizer._configure()
    return
  

  def test_constructor(self):
    """Test constructor.
    """
    output = OutputManager()
    output.inventory.writer._configure()
    output._configure()
    return


  def test_preinitialize(self):
    """Test preinitialize().
    """
    dataProvider = TestProvider()
    output = OutputManager()
    output.preinitialize(dataProvider)
    
    self.failIf(output.dataProvider is None)
    return


  def test_verifyConfiguration(self):
    """Test verifyConfiguration().
    """
    dataProvider = TestProvider()
    output = OutputManager()
    output.preinitialize(dataProvider)

    output.vertexInfoFields = ["vertex info"]
    output.vertexDataFields = ["vertex data 2"]
    output.cellInfoFields = []
    output.cellDataFields = ["cell data"]
    output.verifyConfiguration(dataProvider.getDataMesh())
    return
  
  
  def test_initialize(self):
    """Test initialize().
    """
    # No quadrature
    output = OutputManager()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    # With quadrature
    from pylith.feassemble.FIATLagrange import FIATLagrange
    from pylith.feassemble.Quadrature import Quadrature
    cell = FIATLagrange()
    cell.inventory.dimension = 2
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    
    output = OutputManager()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer, quadrature)
    return


  def test_openclose(self):
    """Test open() and close().
    """
    output = OutputManager()
    output.inventory.writer.inventory.filename = "output.vtk"
    output.inventory.writer._configure()
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    output.open(totalTime=5.0, numTimeSteps=2)
    output.close()
    return


  def test_writeInfo(self):
    """Test writeInfo().
    """
    output = OutputManager()
    output.inventory.writer.inventory.filename = "output.vtk"
    output.inventory.writer._configure()
    output.inventory.vertexInfoFields = ["vertex info"]
    output.inventory.cellInfoFields = ["cell info"]
    output._configure()
    
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)
    
    output.open(totalTime=5.0, numTimeSteps=2)
    output.writeInfo()
    output.close()
    return


  def test_writeData(self):
    """Test writeData().
    """
    output = OutputManager()
    output.inventory.writer.inventory.filename = "output.vtk"
    output.inventory.writer.inventory.timeFormat = "%3.1f"
    output.inventory.writer._configure()
    output.inventory.vertexDataFields = ["vertex data 2",
                                         "vertex data 1"]
    output.inventory.cellDataFields = ["cell data"]
    output._configure()
    
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    output.open(totalTime=5.0, numTimeSteps=2)
    output.writeData(2.0, dataProvider.fields)
    output.close()
    return


  def test_checkWrite(self):
    """Test _checkWrite().
    """
    dataProvider = TestProvider()

    # Default values should be true
    output = OutputManager()
    output.inventory.writer._configure()
    output._configure()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)
    self.assertEqual(True, output._checkWrite(0.0))
    self.assertEqual(True, output._checkWrite(3.234e+8))

    # Check writing based on time
    output = OutputManager()
    output.inventory.writer._configure()
    output._configure()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    output.inventory.outputFreq = "time_step"
    t = 0.0
    dt = output.dtN
    self.assertEqual(True, output._checkWrite(t))
    self.assertEqual(False, output._checkWrite(t))
    self.assertEqual(False, output._checkWrite(t + 0.8*dt))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    t = 2*dt
    self.assertEqual(True, output._checkWrite(t))
    
    # Check writing based on number of steps
    output = OutputManager()
    output.inventory.writer._configure()
    output.inventory.outputFreq = "skip"
    output.inventory.skip = 1
    output._configure()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)
    t = 0.0
    dt = 1.0
    self.assertEqual(True, output._checkWrite(t))
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    output.inventory.skip = 2
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    
    return


  def test_factory(self):
    """Test factory method.
    """
    from pylith.meshio.OutputManager import output_manager
    o = output_manager()
    return


# End of file 
