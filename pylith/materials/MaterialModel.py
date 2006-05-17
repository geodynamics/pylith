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

## @file pylith/materials/MaterialModel.py
## @brief Python abstract base class for material constitutive relation.

from pyre.components.Component import Component

# MaterialModel class
class MaterialModel(Component):
  """Python abstract base class for material constitutive relation."""

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def elasticityConsts(self):
    raise NotImplementedError, \
          "MaterialModel::elasticityConsts() not implemented."
    return

  def __init__(self, name="materialmodel"):
    """Constructor."""
    Component.__init__(self, name, facility="materialmodel")
    self.handle = None
    self.queryVals = []
    return
  
 # version
__id__ = "$Id$"

# End of file 
