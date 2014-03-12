#!/usr/bin/env python
#
# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#

## @file tests/3d/tet4/TestShearDispNoSlipRefine.py
##
## @brief Test suite for testing pylith with 3-D shear motion with no
## fault slip and mesh refinement.

import numpy
from TestTet4 import TestTet4
from sheardisp_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class ShearApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="sheardispnosliprefine")
    return


# Helper function to run PyLith
def run_pylith():
  """
  Run pylith.
  """
  if not "done" in dir(run_pylith):
    # Generate spatial databases
    from sheardisp_gendb import GenerateDB
    db = GenerateDB()
    db.run()

    # Run PyLith
    app = ShearApp()
    run_pylith.done = True # Put before run() so only called once
    app.run()
  return


class TestShearDispNoSlipRefine(TestTet4):
  """
  Test suite for testing pylith with 2-D shear extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTet4.setUp(self)
    self.mesh = {'ncells-elastic': 7068*8,
                 'ncells-viscoelastic': 11131*8,
                 'ncorners': 4,
                 'nvertices': 26947,
                 'spaceDim': 3,
                 'tensorSize': 6}
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 0
    self.faultMesh = {'nvertices': 51,
                      'spaceDim': 3,
                      'ncells': 30*4,
                      'ncorners': 3}
    run_pylith()
    self.outputRoot = "sheardispnosliprefine"

    self.soln = AnalyticalSoln()
    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices)


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    ncells = self.mesh['ncells']
    pts = numpy.zeros( (ncells, 3), dtype=numpy.float64)
    if name == "total_strain":
      stateVar = self.soln.strain(pts)
    elif name == "stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


  def calcFaultField(self, name, vertices):
    """
    Calculate fault info.
    """

    normalDir = (-1.0, 0.0, 0.0)
    finalSlip = 0.0
    slipTime = 0.0

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]
      field[0,:,2] = normalDir[2]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, 3), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestShearDispNoSlipRefine import TestShearDispNoSlipRefine as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
