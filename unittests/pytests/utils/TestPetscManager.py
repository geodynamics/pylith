#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/utils/TestPescManager.py

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
    manager = PetscManager()
    manager.initialize()
    manager.finalize()
    return


# End of file 
