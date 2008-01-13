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

    from pylith.meshio.SolutionIOVTK import SolutionIOVTK
    io = SolutionIOVTK()
    io._configure()
    io.filename = 'mesh.vtk'
    from spatialdata.geocoords.CSCart import CSCart
    io.coordsys = CSCart()
    mesh.coordsys = CSCart()
    io.open(mesh)
    from pyre.units.time import s
    t = 0.0*s
    io.openTimeStep(0.0, 0, mesh)
    io.closeTimeStep()
    io.close()

    from pylith.topology.Distributor import Distributor
    distributor = Distributor()
    distributor.partitioner = "chaco"
    newMesh = distributor.distribute(mesh)
    io.filename = 'newMesh.vtk'
    io.open(newMesh)
    io.openTimeStep(0.0, 0, mesh)
    io.closeTimeStep()
    io.close()
    return


# End of file 
