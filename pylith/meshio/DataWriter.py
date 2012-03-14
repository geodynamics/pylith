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

## @file pyre/meshio/DataWriter.py
##
## @brief Python abstract base class for writing finite-element data.
##
## Factory: output_data_writer

from pylith.utils.PetscComponent import PetscComponent

# DataWriter class
class DataWriter(PetscComponent):
  """
  Python abstract base class for writing finite-element data.
  """

  # INVENTORY //////////////////////////////////////////////////////////
  
  # None

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="datawriter"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="datawriter")
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self, normalizer):
    """
    Initialize writer.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    self.deallocate()
    return


# End of file
