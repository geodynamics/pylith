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


# End of file
