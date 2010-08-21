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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/utils/NullComponent.py
##
## @brief Python NullComponent object that is an empty component.

from PetscComponent import PetscComponent

# NullComponent class
class NullComponent(PetscComponent):
  """
  Python NullComponent object that is an empty component.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name="nullcomponent", facility="nullcomponent")
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    return
    

# End of file 
