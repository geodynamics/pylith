#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/meshio/TestSingleOutput.py

# @brief Unit testing of Homogenous object.

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
