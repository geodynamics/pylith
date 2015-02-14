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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/topology/TestRefineUniform.py

## @brief Unit testing of Python RefineUniform object.

import unittest

from pylith.topology.RefineUniform import RefineUniform

# ----------------------------------------------------------------------
class TestRefineUniform(unittest.TestCase):
  """
  Unit testing of Python RefineUniform object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    io = RefineUniform()
    return


  def test_refineTet4NoFault(self):
    """
    Test refine().
    """
    filenameIn = "data/twotet4.mesh"
    filenameOut = "data/twotet4_test.mesh"
    filenameOutE = "data/twotet4_nofault_refined2.mesh"

    self._runTest(filenameIn, filenameOut, filenameOutE)
    return


  def test_refineTet4Fault(self):
    """
    Test refine().
    """
    filenameIn = "data/twotet4.mesh"
    filenameOut = "data/twotet4_test.mesh"
    filenameOutE = "data/twotet4_fault_refined2.mesh"

    self._runTest(filenameIn, filenameOut, filenameOutE, "fault")
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.topology.RefineUniform import mesh_refiner
    refiner = mesh_refiner()
    return


  def _runTest(self, filenameIn, filenameOut, filenameOutE, faultGroup=None):

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs._configure()

    from pylith.meshio.MeshIOAscii import MeshIOAscii
    io = MeshIOAscii()
    io.inventory.filename = filenameIn
    io.inventory.coordsys = cs
    io._configure()
    
    mesh = io.read(debug=False, interpolate=True)

    if not faultGroup is None:
      from pylith.faults.FaultCohesiveKin import FaultCohesiveKin
      fault = FaultCohesiveKin()
      fault.inventory.matId = 10
      fault.inventory.faultLabel = faultGroup
      fault._configure()

      nvertices = fault.numVerticesNoMesh(mesh)
      firstFaultVertex = 0
      firstLagrangeVertex = nvertices
      firstFaultCell      = 2*nvertices
      fault.adjustTopology(mesh, 
                           firstFaultVertex, 
                           firstLagrangeVertex,
                           firstFaultCell)

    from pylith.topology.RefineUniform import RefineUniform
    refiner = RefineUniform()
    meshRefined = refiner.refine(mesh)

    return



# End of file 
