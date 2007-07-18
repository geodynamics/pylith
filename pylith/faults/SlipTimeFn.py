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

## @file pylith/faults/SlipTimeFn.py
##

## @brief Python abstract base class for kinematic slip time function.
##
## Factory: slip_time_fn

from pyre.components.Component import Component

# SlipTimeFn class
class SlipTimeFn(Component):
  """
  Python abstract base class for kinematic slip time function.

  Factory: slip_time_fn
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="sliptimefn"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="sliptimefn")
    self.cppHandle = None
    return


  def preinitialize(self):
    """
    Do pre-initialization setup.
    """
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    return


  def initialize(self):
    """
    Initialize.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    return

  
# End of file 
