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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/3dnew/hex8/TestHex8.py
##
## @brief Generic tests for problems using 3-D mesh.

import unittest
import numpy

class TestHex8(unittest.TestCase):
  """
  Generic tests for problems using 3-D mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.mesh = {'ncells_elastic': 128,
                 'ncells_viscoelastic': 256,
                 'ncorners': 8,
                 'nvertices': 567,
                 'spaceDim': 3,
                 'tensorSize': 6}
    return


  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if self.reader is None:
      return

    from pylith.tests.PhysicalProperties import check_properties
    from rigidbody_soln import p_mu,p_lambda,p_density

    self.mesh['ncells'] = self.mesh['ncells_elastic']
    ncells= self.mesh['ncells']
    filename = "%s-elastic_info.vtk" % self.outputRoot
    propMu =  p_mu*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propLambda = p_lambda*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propDensity = p_density*numpy.ones( (ncells, 2), dtype=numpy.float64)
    properties = {'mu': propMu,
                  'lambda': propLambda,
                  'density': propDensity}
    check_properties(self, filename, self.mesh, properties)

    self.mesh['ncells'] = self.mesh['ncells_viscoelastic']
    ncells= self.mesh['ncells']
    filename = "%s-viscoelastic_info.vtk" % self.outputRoot
    propMu =  p_mu*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propLambda = p_lambda*numpy.ones( (ncells, 1), dtype=numpy.float64)
    propDensity = p_density*numpy.ones( (ncells, 2), dtype=numpy.float64)
    properties = {'mu': propMu,
                  'lambda': propLambda,
                  'density': propDensity}
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

    from pylith.tests.StateVariables import check_state_variables
    stateVars = ["total_strain", "stress"]

    filename = "%s-elastic_t0000000.vtk" % self.outputRoot
    self.mesh['ncells'] = self.mesh['ncells_elastic']
    check_state_variables(self, filename, self.mesh, stateVars)

    filename = "%s-viscoelastic_t0000000.vtk" % self.outputRoot
    self.mesh['ncells'] = self.mesh['ncells_viscoelastic']
    check_state_variables(self, filename, self.mesh, stateVars)

    return


# End of file
