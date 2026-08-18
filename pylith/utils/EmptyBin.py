# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.petsc.Component import Component


class EmptyBin(Component):
  """
  Empty container for a collection of objects.
  """

  def __init__(self, name="emptybin"):
    """Constructor.
    """
    Component.__init__(self, name, facility="empty_bin")


# End of file
