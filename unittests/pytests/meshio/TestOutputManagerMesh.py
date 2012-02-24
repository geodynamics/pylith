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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/meshio/TestOutputManagerMesh.py

## @brief Unit testing of Python OutputManagerMesh object.

import unittest

from pylith.meshio.OutputManagerMesh import OutputManagerMesh

# ----------------------------------------------------------------------
class TestProvider(object):

  def __init__(self):
    """
    Constructor.
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

    from pylith.topology.Fields import MeshFields
    fields = MeshFields(mesh)
    
    self.mesh = mesh
    self.fields = fields
    return


  def getDataMesh(self):
    """
    Get mesh.
    """
    return (self.mesh, None, None)


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    from pylith.field.field.FieldBase import SCALAR,VECTOR,OTHER
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
    """
    Get cell field.
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
class TestOutputManagerMesh(unittest.TestCase):
  """
  Unit testing of Python OutputManagerMesh object.
  """

  def setUp(self):
    from spatialdata.units.Nondimensional import Nondimensional
    self.normalizer = Nondimensional()
    self.normalizer._configure()
    return
  

  def test_constructor(self):
    """
    Test constructor.
    """
    output = OutputManagerMesh()
    output.inventory.writer._configure()
    output._configure()
    return


  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    dataProvider = TestProvider()
    output = OutputManagerMesh()
    output.preinitialize(dataProvider)
    
    self.assertEqual(dataProvider, output.dataProvider)
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().
    """
    dataProvider = TestProvider()
    output = OutputManagerMesh()
    output.preinitialize(dataProvider)

    output.vertexInfoFields = ["vertex info"]
    output.vertexDataFields = ["vertex data 2"]
    output.cellInfoFields = []
    output.cellDataFields = ["cell data"]
    output.verifyConfiguration(dataProvider.getDataMesh())
    return
  
  
  def test_initialize(self):
    """
    Test initialize().
    """
    # No quadrature
    output = OutputManagerMesh()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    # With quadrature
    from pylith.feassemble.FIATLagrange import FIATLagrange
    from pylith.feassemble.Quadrature import MeshQuadrature
    cell = FIATLagrange()
    cell.inventory.dimension = 3
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    quadrature = MeshQuadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()
    
    output = OutputManagerMesh()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer, quadrature)
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    output = OutputManagerMesh()
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
    """
    Test writeInfo().
    """
    output = OutputManagerMesh()
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
    """
    Test writeData().
    """
    output = OutputManagerMesh()
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


  def test_estimateNumSteps(self):
    """
    Test _estimateNumSteps().
    """
    from pyre.units.time import second

    output = OutputManagerMesh()
    output.inventory.outputFreq = "skip"
    output.inventory.skip = 2
    output._configure()

    numTimeSteps = 0
    totalTime = 1.0*second
    self.assertEqual(0, output._estimateNumSteps(totalTime, numTimeSteps))

    numTimeSteps = 4
    totalTime = 1.0*second
    self.assertEqual(2, output._estimateNumSteps(totalTime, numTimeSteps))

    numTimeSteps = 2
    totalTime = 1.0*second
    self.assertEqual(1, output._estimateNumSteps(totalTime, numTimeSteps))

    output = OutputManagerMesh()
    output.inventory.outputFreq = "time_step"
    output.inventory.dt = 2.0*second
    output._configure()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)

    numTimeSteps = 0
    totalTime = 1.0*second
    self.assertEqual(0, output._estimateNumSteps(totalTime, numTimeSteps))

    numTimeSteps = 4
    totalTime = 1.0*second
    self.assertEqual(1, output._estimateNumSteps(totalTime, numTimeSteps))

    numTimeSteps = 2
    totalTime = 2.0*second
    self.assertEqual(2, output._estimateNumSteps(totalTime, numTimeSteps))
    return


  def test_checkWrite(self):
    """
    Test _checkWrite().
    """
    dataProvider = TestProvider()

    # Default values should be true
    output = OutputManagerMesh()
    output.inventory.writer._configure()
    output._configure()
    output.preinitialize(dataProvider)
    output.initialize(self.normalizer)
    self.assertEqual(True, output._checkWrite(0.0))
    self.assertEqual(True, output._checkWrite(3.234e+8))

    # Check writing based on time
    output = OutputManagerMesh()
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
    output = OutputManagerMesh()
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
    """
    Test factory method.
    """
    from pylith.meshio.OutputManagerMesh import output_manager
    o = output_manager()
    return


# End of file 
