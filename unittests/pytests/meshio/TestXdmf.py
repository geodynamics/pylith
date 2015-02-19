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

## @file unittests/pytests/meshio/TestXdmf.py

## @brief Unit testing of Python Xdmf object.

import unittest

from pylith.meshio.Xdmf import Xdmf

# ----------------------------------------------------------------------
class TestXdmf(unittest.TestCase):
  """
  Unit testing of Python OutputManagerMesh object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    writer = Xdmf()
    writer._configure()
    return


# End of file 
