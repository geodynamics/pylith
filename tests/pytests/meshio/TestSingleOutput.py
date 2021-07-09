#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#

## @file tests/pytests/meshio/TestSingleOutput.py

## @brief Unit testing of Homogenous object.

import unittest

# ----------------------------------------------------------------------
class TestSingleOutput(unittest.TestCase):
  """Unit testing of SingleOutput object.
  """

  def test_constructor(self):
    """Test constructor.
    """
    from pylith.meshio.SingleOutput import SingleOutput
    outputs = SingleOutput()
    return


  def test_configure(self):
    """Test _configure().
    """
    from pylith.meshio.SingleOutput import SingleOutput
    outputs = SingleOutput()
    from pylith.meshio.OutputSoln import OutputSoln
    outputs.inventory.output = OutputSoln()
    outputs._configure()
    self.assertEqual(1, len(outputs.components()))
    return


# End of file 
