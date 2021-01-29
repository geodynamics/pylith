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

## @file tests/pytests/meshio/TestOutputSolnPoints.py

## @brief Unit testing of Python OutputSolnPoints object.

import unittest

from pylith.meshio.OutputSolnPoints import OutputSolnPoints

# ----------------------------------------------------------------------
class TestOutputSolnPoints(unittest.TestCase):
  """
  Unit testing of Python OutputSolnPoints object.
  """

  def setUp(self):
    from pylith.meshio.MeshIOAscii import MeshIOAscii
    iohandler = MeshIOAscii()
    filename = "data/twohex8.txt"
    
    from spatialdata.units.Nondimensional import Nondimensional
    normalizer = Nondimensional()
    normalizer._configure()

    from spatialdata.geocoords.CSCart import CSCart
    iohandler.inventory.filename = filename
    iohandler.inventory.coordsys = CSCart()
    iohandler._configure()
    mesh = iohandler.read(debug=False, interpolate=False)

    from pylith.topology.SolutionFields import SolutionFields
    fields = SolutionFields(mesh)

    name = "disp(t)"
    fields.add(name, "displacement")
    fields.solutionName(name)
    field = fields.get(name)
    field.newSection(field.VERTICES_FIELD, mesh.dimension())
    field.allocate()

    self.mesh = mesh
    self.fields = fields
    self.normalizer = normalizer
    return


  def test_constructor(self):
    """
    Test constructor.
    """
    output = OutputSolnPoints()
    output.inventory.writer._configure()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output._configure()
    return


  def test_preinitialize(self):
    """
    Test preinitialize().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output._configure()
    output.preinitialize()
    
    self.failIf(output.dataProvider is None)
    return


  def test_verifyConfiguration(self):
    """
    Test verifyConfiguration().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output._configure()
    output.preinitialize()

    output.vertexDataFields = ["displacement"]
    output.verifyConfiguration(self.mesh)
    return
  
  
  def test_initialize2(self):
    """
    Test initialize().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/point.txt"
    output.inventory.reader._configure()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()

    output.preinitialize()
    output.initialize(self.mesh, self.normalizer)
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()

    output.preinitialize()
    output.initialize(self.mesh, self.normalizer)
    return


  def test_openclose(self):
    """
    Test open() and close().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output.inventory.writer.inventory.filename = "test.vtk"
    output.inventory.writer._configure()
    output._configure()

    output.preinitialize()
    output.initialize(self.mesh, self.normalizer)

    from pythia.pyre.units.time import s
    output.open(totalTime=5.0*s, numTimeSteps=2)
    output.close()
    return


  def test_writeInfo(self):
    """
    Test writeInfo().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output.inventory.writer.inventory.filename = "output_sub.vtk"
    output.inventory.writer._configure()
    output._configure()

    output.preinitialize()
    output.initialize(self.mesh, self.normalizer)
    
    output.open(totalTime=5.0, numTimeSteps=2)
    output.writeInfo()
    output.close()
    return


  def test_writeData(self):
    """
    Test writeData().
    """
    output = OutputSolnPoints()
    output.inventory.reader.inventory.filename = "data/points.txt"
    output.inventory.reader._configure()
    output.inventory.writer.inventory.filename = "output_sub.vtk"
    output.inventory.writer.inventory.timeFormat = "%3.1f"
    output.inventory.writer._configure()
    output.inventory.vertexDataFields = ["displacement"]
    output._configure()

    output.preinitialize()
    output.initialize(self.mesh, self.normalizer)

    output.open(totalTime=5.0, numTimeSteps=2)
    output.writeData(2.0, self.fields)
    output.close()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.meshio.OutputSolnPoints import output_manager
    o = output_manager()
    return


# End of file 
