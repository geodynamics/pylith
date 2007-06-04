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
    io.filename = 'mesh.vtk'
    from spatialdata.geocoords.CSCart import CSCart
    io.coordsys = CSCart()
    mesh.coordsys = CSCart()
    io.open(mesh)
    io.writeTopology()
    io.close()

    #from pylith.topology.Distributor import Distributor
    #distributor = Distributor()
    #newMesh = distributor.distribute(mesh)
    #io.filename = 'newMesh.vtk'
    #io.open(newMesh)
    #io.writeTopology(newMesh)
    #io.close()
    return


# End of file 
