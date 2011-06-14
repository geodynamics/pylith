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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/tri3/TestTri3.py
##
## @brief Generic tests for problems using 2-D mesh.

import unittest
import numpy

class TestTri3(unittest.TestCase):
  """
  Generic tests for problems using 2-D mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.mesh = {'ncells': 124,
                 'ncorners': 3,
                 'nvertices': 79,
                 'spaceDim': 3,
                 'tensorSize': 3}
    return


  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if self.reader is None:
      return

    ncells= self.mesh['ncells']

    filename = "%s-elastic_info.vtk" % self.outputRoot
    from axialdisp_soln import p_mu,p_lambda,p_density

    propMu =  p_mu*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propLambda = p_lambda*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propDensity = p_density*numpy.ones( (ncells, 2), dtype=numpy.float64)

    properties = {'mu': propMu,
                  'lambda': propLambda,
                  'density': propDensity}

    from pylith.tests.PhysicalProperties import check_properties
    check_properties(self, filename, self.mesh, properties)

    return


  def test_soln(self):
    """
    Check solution (displacement) field.
    """
    if self.reader is None:
      return

    filename = "%s_t0000000.vtk" % self.outputRoot
    from pylith.tests.Solution import check_displacements
    check_displacements(self, filename, self.mesh)

    return


  def test_elastic_statevars(self):
    """
    Check elastic state variables.
    """
    if self.reader is None:
      return

    filename = "%s-elastic_t0000000.vtk" % self.outputRoot

    from pylith.tests.StateVariables import check_state_variables
    stateVars = ["total_strain", "stress"]
    check_state_variables(self, filename, self.mesh, stateVars)

    return


# End of file
