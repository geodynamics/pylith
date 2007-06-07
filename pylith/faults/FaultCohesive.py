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

## @file pylith/faults/FaultCohesive.py
##

## @brief Python abstract base class for a fault surface implemented
## with cohesive elements.
##
## Factory: fault

from Fault import Fault

# FaultCohesive class
class FaultCohesive(Fault):
  """
  Python abstract base class for a fault surface implemeted with
  cohesive elements.

  Factory: fault
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesive"):
    """
    Constructor.
    """
    Fault.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Fault._configure(self)
    return

  
# End of file 
