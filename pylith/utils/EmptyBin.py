# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#
# @file pylith/utils/EmptyBin.py
#
# @brief Python container for a collection of objects.
#
# Factory: object_bin

from pylith.utils.PetscComponent import PetscComponent


class EmptyBin(PetscComponent):
  """Python container for a collection of objects.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="emptybin"):
    """Constructor.
    """
    PetscComponent.__init__(self, name, facility="empty_bin")
    return


# End of file
