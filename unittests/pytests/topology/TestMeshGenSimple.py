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
    io.initialize()
    io.open(mesh)
    from pyre.units.time import s
    t = 0.0*s
    io.openTimeStep(t, 0)
    io.closeTimeStep()
    io.close()

    from pylith.topology.Distributor import Distributor
    distributor = Distributor()
    distributor.partitioner = "chaco"
    newMesh = distributor.distribute(mesh)
    io.writer.filename = 'newMesh.vtk'
    io.open(newMesh)
    io.openTimeStep(t, 0)
    io.closeTimeStep()
    io.close()
    return


# End of file 
