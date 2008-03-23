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

## @file unittests/pytests/meshio/TestOutputManager.py

## @brief Unit testing of Python OutputManager object.

import unittest

from pylith.meshio.OutputManager import OutputManager

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

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    iohandler = MeshIOAscii()
    filename = "data/twohex8.txt"
    
    from spatialdata.geocoords.CSCart import CSCart
    iohandler.filename = filename
    iohandler.coordsys = CSCart()
    mesh = iohandler.read(debug=False, interpolate=False)

    from pylith.topology.FieldsManager import FieldsManager
    fields = FieldsManager(mesh)
    
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
    if name == "vertex info":
      fieldType = 0
      fiberDim = 1
    elif name == "vertex data 1":
      fieldType = 1
      fiberDim = self.mesh.dimension()
    elif name == "vertex data 2":
      fieldType = 3
      fiberDim = 5
    else:
      raise ValueError("Unknown field '%s'." % name)

    self.fields.addReal(name)
    self.fields.setFiberDimension(name, fiberDim)
    self.fields.allocate(name)
    field = self.fields.getReal(name)
    return (field, fieldType)


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if name == "cell info":
      fieldType = 0
      fiberDim = 1
    elif name == "cell data":
      fieldType = 1
      fiberDim = self.mesh.dimension()
    else:
      raise ValueError("Unknown field '%s'." % name)

    self.fields.addReal(name)
    self.fields.setFiberDimension(name, fiberDim, "cells")
    self.fields.allocate(name)
    field = self.fields.getReal(name)
    return (field, fieldType)


# ----------------------------------------------------------------------
class TestOutputManager(unittest.TestCase):
  """
  Unit testing of Python OutputManager object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    output = OutputManager()
    output._configure()
    output.writer._configure()
    return


  def test_sync(self):
    """
    Test _sync().
    """
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.initialize()
    output._sync()
    self.assertNotEqual(None, output.cppHandle)
    return
  

  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    output = OutputManager()
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    
    self.assertEqual(dataProvider, output.dataProvider)
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().
    """
    output = OutputManager()
    dataProvider = TestProvider()
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
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.filename = "test.vtk"
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize()
    self.assertNotEqual(None, output.cppHandle)

    # With quadrature
    from pylith.feassemble.FIATSimplex import FIATSimplex
    from pylith.feassemble.quadrature.Quadrature1D import Quadrature1D
    cell = FIATSimplex()
    cell.shape = "line"
    cell.degree = 2
    cell.order = 2

    quadrature = Quadrature1D()
    quadrature.cell = cell
    
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.filename = "test.vtk"
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize(quadrature)
    self.assertNotEqual(None, output.cppHandle)
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.filename = "output.vtk"
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize()

    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.close()
    return


  def test_writeInfo(self):
    """
    Test writeInfo().
    """
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.filename = "output.vtk"
    output.vertexInfoFields = ["vertex info"]
    output.cellInfoFields = ["cell info"]
    
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize()
    
    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.writeInfo()
    output.close()
    return


  def test_writeData(self):
    """
    Test writeData().
    """
    output = OutputManager()
    output._configure()
    output.writer._configure()
    output.writer.filename = "output.vtk"
    output.writer.timeFormat = "%3.1f"
    output.vertexDataFields = ["vertex data 2",
                               "vertex data 1"]
    output.cellDataFields = ["cell data"]
    
    dataProvider = TestProvider()
    output.preinitialize(dataProvider)
    output.initialize()

    from pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.writeData(2.0*s, dataProvider.fields)
    output.close()
    return


  def test_checkWrite(self):
    """
    Test _checkWrite().
    """
    from pyre.units.time import s

    # Default values should be true
    output = OutputManager()
    output._configure()
    self.assertEqual(True, output._checkWrite(0.0*s))
    self.assertEqual(True, output._checkWrite(3.234e+8*s))

    # Check writing based on time
    output = OutputManager()
    output._configure()
    output.outputFreq = "time_step"
    t = 0.0*s
    dt = output.dt
    self.assertEqual(True, output._checkWrite(t))
    self.assertEqual(False, output._checkWrite(t))
    self.assertEqual(False, output._checkWrite(t + 0.8*dt))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    t = 2*dt
    self.assertEqual(True, output._checkWrite(t))
    
    # Check writing based on number of steps
    output = OutputManager()
    output._configure()
    output.outputFreq = "skip"
    output.skip = 1
    t = 0.0*s
    dt = 1.0*s
    self.assertEqual(True, output._checkWrite(t))
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    output.skip = 2
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(False, output._checkWrite(t))
    t += dt
    self.assertEqual(True, output._checkWrite(t))
    
    return


# End of file 
