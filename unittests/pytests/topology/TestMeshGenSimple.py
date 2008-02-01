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

## @file unittests/pytests/topology/TestMeshGenSimple.py

## @brief Unit testing of MeshGenSimple object.

import unittest

from pylith.topology.MeshGenSimple import MeshGenSimple

# ----------------------------------------------------------------------
class TestMeshGenSimple(unittest.TestCase):
  """
  Unit testing of MeshGenSimple object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    generator = MeshGenSimple()
    return
  

  def test_distribute(self):
    """
    Test distribute()
    """
    generator = MeshGenSimple()
    generator.interpolate = True
    generator.setBoundary(generator.createCubeBoundary())
    mesh = generator.create()

    from pylith.meshio.OutputManager import OutputManager
    io = OutputManager()
    io._configure()
    io.writer._configure()
    io.writer.filename = 'mesh.vtk'
    from spatialdata.geocoords.CSCart import CSCart
    io.coordsys = CSCart()
    mesh.coordsys = CSCart()
    self.mesh = mesh
    
    from pyre.units.time import s
    t = 0.0*s
    io.preinitialize(self)
    io.initialize()
    io.writeData(t)

    from pylith.topology.Distributor import Distributor
    distributor = Distributor()
    distributor.partitioner = "chaco"
    newMesh = distributor.distribute(mesh)
    self.mesh = newMesh
    io.writer.filename = 'newMesh.vtk'
    io.writeData(t)
    return


  def getDataMesh(self):
    """
    Get mesh.
    """
    return (self.mesh, None, None)
  

# End of file 
