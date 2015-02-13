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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/topology/SolutionFields.py
##
## @brief Python object for managing fields associated with problem
## solution.

from topology import SolutionFields as ModuleSolutionFields

# ----------------------------------------------------------------------
# SolutionFields class
class SolutionFields(ModuleSolutionFields):
  """
  Python object for managing fields associated with problem solution.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    ModuleSolutionFields.__init__(self, mesh)
    return
    

  def cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    self.deallocate()
    return
    

# End of file
