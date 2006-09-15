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

## @file pylith/topology/MeshGenerator.py
## @brief Python abstract base class for mesh generator.

from pyre.components.Component import Component

# MeshGenerator class
class MeshGenerator(Component):
  """Python abstract base class for mesh generator."""

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def create(self):
    raise NotImplementedError, "MeshGenerator::create() not implemented."
    return


  def __init__(self, name="meshgenerator"):
    """Constructor."""
    Component.__init__(self, name, facility="meshgenerator")
    return


# version
__id__ = "$Id$"

# End of file 
