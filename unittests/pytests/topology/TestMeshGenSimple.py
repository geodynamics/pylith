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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
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
    io.preinitialize(self)
    io.initialize()
    io.writeInfo()

    from pylith.topology.Distributor import Distributor
    distributor = Distributor()
    distributor.partitioner = "chaco"
    newMesh = distributor.distribute(mesh)
    self.mesh = newMesh
    io.writer.filename = 'newMesh.vtk'
    io.writeInfo()
    return


  def getDataMesh(self):
    """
    Get mesh.
    """
    return (self.mesh, None, None)
  

# End of file 
