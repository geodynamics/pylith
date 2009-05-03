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

## @file pylith/problems/TimeStepAdapt.py
##
## @brief Python class for marching forward in time with irregular
## time steps that adapt to the solution.
##
## Factory: time_step

from TimeStep import TimeStep

from pylith.utils.profiling import resourceUsageString

# TimeStepAdapt class
class TimeStepAdapt(TimeStep):
  """
  Python abstract base class for marching format in time with
  irregular time steps that adapt to the solution.

  Factory: time_step.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(TimeStep.Inventory):
    """
    Python object for managing TimeStepAdapt facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing TimeStepAdapt facilities and properties.
    ##
    ## \b Properties
    ## @li \b total_time Time duration for simulation.
    ## @li \b max_dt Maximum time step.
    ## @li \b adapt_skip Number of time steps to skip between adjusting value.
    ## @li \b stability_factor "Safety factor" for stable time step.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    maxDt = pyre.inventory.dimensional("max_dt", default=1.0*second,
                                    validator=pyre.inventory.greater(0.0*second))
    maxDt.meta['tip'] = "Maximum time step permitted."

    adaptSkip = pyre.inventory.int("adapt_skip", default=10,
                                   validator=pyre.inventory.greaterEqual(0))
    adaptSkip.meta['tip'] = "Number of time steps to skip between " \
        "adjusting value."

    stabilityFactor = pyre.inventory.float("stability_factor", default=1.2,
                                    validator=pyre.inventory.greater(0.0))
    stabilityFactor.meta['tip'] = "'Safety factor' for stable time step."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestepadapt"):
    """
    Constructor.
    """
    TimeStep.__init__(self, name)
    self._loggingPrefix = "DtAd "
    self.skipped = 0
    return


  def initialize(self, normalizer):
    """
    Initialize time step algorithm.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    TimeStep.initialize(self, normalizer)

    # Nondimensionalize time scales
    timeScale = normalizer.timeScale()
    self.maxDt = normalizer.nondimensionalize(self.maxDt, timeScale)

    self._logger.eventEnd(logEvent)
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    # Guess using maximum time step
    nsteps = int(1.0 + self.totalTime / self.maxDt)
    return nsteps


  def timeStep(self, mesh, integrators):
    """
    Adjust stable time step for advancing forward in time.
    """
    dtStable = 1.0e+30
    for integrator in integrators:
      dt = integrator.stableTimeStep(mesh)
      if dt < dtStable:
        dtStable = dt
    
    if self.skipped < self.adaptSkip and \
          self.dt != 0.0 and \
          self.dt < dtStable:
      self.skipped += 1
    else:
      self.dt = min(dtStable/self.stabilityFactor, self.maxDt)
      self.skipped = 0
    return self.dt

  
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    TimeStep._configure(self)
    self.totalTime = self.inventory.totalTime
    self.maxDt = self.inventory.maxDt
    self.adaptSkip = self.inventory.adaptSkip
    self.stabilityFactor = self.inventory.stabilityFactor
    self.dt = self.maxDt
    return


# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStepAdapt.
  """
  return TimeStepAdapt()


# End of file 
