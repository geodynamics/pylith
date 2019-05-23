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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
##
## @file tests_auto/linearelasticity/nofaults/TestTri.py
##
## @brief Base class for tests with 2-D simplex mesh.

import unittest
import numpy

from pylith.tests import has_h5py

class TestTri3(unittest.TestCase):
  """
  Generic tests for problems using 2-D mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.spaceDim = 2
    self.tensorSize = 3
    self.domain = {
      "ncells": 142,
      "ncorners": 3,
      "nvertices": 88,
    }
    self.materials = {
      "elastic_xpos": {
        "ncells": 0,
        "ncorners": 3,
        "nvertices": 0,
      },
      "elastic_xneg": {
        "ncells": 0,
        "ncorners": 3,
        "nvertices": 0,
      }
    }

    if has_h5py():
      self.checkResults = True
    else:
      self.checkResults = False
    return


  def test_solution_domain(self):
    """
    Check solution (displacement) field.
    """
    if not self.checkResults:
      return

    from pylith.tests.Solution import check_vertex_field
    filename = "%s.h5" % self.outputRoot
    check_vertex_field(self, "displacement", filename, self.domain)
    return


  def test_material_info(self):
    """
    Check elastic info.
    """
    if not self.checkResults:
      return

    from pylith.tests.PhysicalProperties import check_properties
    from axialdisp_soln import p_mu,p_lambda,p_density

    for name, matinfo in self.materials.items():
      filename = "{0}-{1}_info.h5".format(self.outputRoot, name)
      ncells= matinfo["ncells"]
      cell_fields = {
        "mu": p_mu*numpy.ones( (1, ncells, 1), dtype=numpy.float64),
        "lambda": p_lambda*numpy.ones( (1, ncells, 1), dtype=numpy.float64),
        "density": p_density*numpy.ones( (1, ncells, 1), dtype=numpy.float64),
      }
      check_cell_info_fields(self, filename, matinfo, properties, self.spaceDim)
    return


  def test_material_solution(self):
    """
    Check material solution.
    """
    if not self.checkResults:
      return

    from pylith.tests.StateVariables import check_state_variables

    filename = "%s-elastic.h5" % self.outputRoot
     = ["total_strain", "stress", "cauchy_stress"]
    check_vertex_fields(self, filename, matinfo, )
    return


# End of file
