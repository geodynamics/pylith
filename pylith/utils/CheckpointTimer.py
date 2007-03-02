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

## @file pylith/utils/CheckpointTimer.py
##
## @brief Python CheckpointTimer object for managing checkpointing.
##
## USAGE:
##
## @li Call initialize with argument 'toplevel' set to top-level
## object that contains a checkpoint() method.
##
## @li Call update() every time step to checkpoint at desired frequency.
##
## Factory: checkpointer.

from pyre.components.Component import Component

# CheckpointTimer class
class CheckpointTimer(Component):
  """
  Python CheckpointTimer object for managing checkpointing.

  USAGE:

  (1) Call initialize with argument 'toplevel' set to top-level object
  that contains a checkpoint() method.

  (2) Call update() every time step to checkpoint at desired frequency.

  Factory: checkpointer.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing CheckpointTimer facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing CheckpointTimer facilities and properties.
    ##
    ## \b Properties
    ## @li dt Simulation time between checkpoints.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    dt = pyre.inventory.dimensional("dt", default=9.9e+99*second,
                          validator=pyre.inventory.greater(0.0*second))
    dt.meta['tip'] = "Simulation time between checkpoints."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="checkpointtimer"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="checkpointer")

    from pyre.units.time import second
    self.t = -8.9e+99*second

    self.toplevel = None
    return


  def initialize(self, toplevel):
    """
    Initialize checkpoint timer.
    """
    self.toplevel = toplevel
    return
  

  def update(self, t):
    """
    CheckpointTimer if necessary.
    """

    if t.value > self.t.value + self.dt.value:
      if app is None:
        raise ValueError, "Atttempting to checkpoint without " \
              "setting toplevel attribute in CheckpointTimer."
      self.toplevel.checkpoint()
      self.t = t
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.dt = self.inventory.dt
    return


# FACTORIES ////////////////////////////////////////////////////////////

def checkpointer():
  """
  Factory associated with CheckpointTimer.
  """
  return CheckpointTimer()


# End of file 
