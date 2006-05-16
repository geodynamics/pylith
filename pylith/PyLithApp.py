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

## @file pylith/PyLithApp.py
## @brief Python PyLith application

from pyre.applications.Script import Script

# PyLithApp class
class PyLithApp(Script):
  """Python PyLithApp application."""
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """Python object for managing PyLithApp facilities and properties."""

    ## @class Inventory
    ## Python object for managing PyLithApp facilities and properties.
    ##
    ## \b Properties
    ## @li totalTime Time duration for simulation
    ##
    ## \b Facilities
    ## @li \b problem Computational problem to solve

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                                  validator=pyre.inventory.greaterEqual(0.0))
    totalTime.meta['tip'] = "Time duration for simulation."

    from pylith.problems.QuasiStatic import QuasiStatic
    problem = pyre.inventory.facility("problem", factory=QuasiStatic)
    problem.meta['tip'] = "Computational problem to solve."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def main(self):
    """Run the application."""

    self.problem.initialize()

    from pyre.units.time import second
    t = 0.0*second
    while t.value < self.totalTime:
      self.problem.prestep()
      dt = self.problem.stableTimeStep()
      self.problem.step(dt)
      self.poststep(t+dt)
      t += dt
    return
  

  def __init__(self, name="pylithapp"):
    """Constructor."""
    Script.__init__(self, name)
    self.totalTime = None
    self.problem = None
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Setup members using inventory."""
    Script._configure(self)
    self.totalTime = self.inventory.totalTime
    self.problem = self.inventory.problem
    return


# version
__id__ = "$Id$"

# End of file 
