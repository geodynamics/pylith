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

## @file unittests/pytests/utils/TestPetscManager.py

## @brief Unit testing of PetscManager object.

import unittest

# ----------------------------------------------------------------------
class TestPetscManager(unittest.TestCase):
  """
  Unit testing of PetscManager object.
  """

  def test_initfinal(self):
    """
    Test initialize/finalize.
    """
    from pylith.utils.PetscManager import PetscManager
    from pylith.utils.petsc import optionsSetValue
    manager = PetscManager()
    manager.options = [("vec_type", "seq")]
    manager.initialize()
    manager.finalize()
    return


# End of file 
