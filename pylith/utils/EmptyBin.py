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

from pylith.utils.PetscComponent import PetscComponent


class EmptyBin(PetscComponent):
  """
  Empty container for a collection of objects.
  """

  def __init__(self, name="emptybin"):
    """Constructor.
    """
    PetscComponent.__init__(self, name, facility="empty_bin")


# End of file
