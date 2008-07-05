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

## @file pylith/problems/TimeStepUniform.py
##
## @brief Python class for marching forward in time with a uniform time step.
##
## Factory: time_step

from TimeStep import TimeStep

from pylith.utils.profiling import resourceUsageString

# TimeStepUniform class
class TimeStepUniform(TimeStep):
  """
  Python abstract base class for marching format in time with a uniform time step.

  Factory: time_step.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(TimeStep.Inventory):
    """
    Python object for managing TimeStepUniform facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing TimeStepUniform facilities and properties.
    ##
    ## \b Properties
    ## @li \b total_time Time duration for simulation.
    ## @li \b default_dt Default time step.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    dt = pyre.inventory.dimensional("dt", default=1.0*second,
                                    validator=pyre.inventory.greater(0.0*second))
    dt.meta['tip'] = "Time step for simulation."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestepuniform"):
    """
    Constructor.
    """
    TimeStep.__init__(self, name)
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    nsteps = 1.0 + int(self.totalTime / self.dt)
    return nsteps


  def timeStep(self, dtStable):
    """
    Adjust stable time step for advancing forward in time.
    """
    return self.dt

  
  def currentStep(self):
    """
    Get current time step size.
    """
    return self.dt
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    TimeStep._configure(self)
    self.totalTime = self.inventory.totalTime
    self.dt = self.inventory.dt
    return


# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStepUniform.
  """
  return TimeStepUniform()


# End of file 
