#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file tests/2d/quad4/TestFrictionQuad4.py
##
## @brief Generic tests for problems using 2-D mesh.

import unittest
import numpy

class TestFrictionQuad4(unittest.TestCase):
  """
  Generic tests for problems using 2-D mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.mesh = {'ncells': 64,
                 'ncorners': 4,
                 'nvertices': 81,
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
    from friction_compression_soln import p_mu,p_lambda,p_density

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
