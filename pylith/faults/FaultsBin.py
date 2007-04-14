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

## @file pylith/faults/FaultsBin.py
##
## @brief Python container for faults.
##
## Factory: faults_bin

from pyre.components.Component import Component

# FaultsBin class
class FaultsBin(Component):
  """
  Python container for faults.

  Factory: faults_bin
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultsbin"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="faults_bin")
    self.faults = []
    return


# End of file 
