#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
  Python object for marching format in time with a uniform time step.

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
    ## @li \b dt Time step for simulation.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    dt = pyre.inventory.dimensional("dt", default=1.0*second,
                                    validator=pyre.inventory.greater(0.0*second))
    dt.meta['tip'] = "Time step for simulation."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestepuniform"):
    """
    Constructor.
    """
    TimeStep.__init__(self, name)
    self._loggingPrefix = "DtUn "
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    nsteps = int(1.0 + self.totalTimeN / self.dtN)
    return nsteps


  def timeStep(self, mesh, integrators):
    """
    Adjust stable time step for advancing forward in time.
    """
    dtStable = self._stableTimeStep(mesh, integrators)
    if dtStable < self.dtN:
      raise RuntimeError("Current nondimensionalized time step of %12.4e "
                         "exceeds the nondimensionalized stable time "
                         "step of %12.4e." % (self.dtN, dtStable))

    return self.dtN

  
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    TimeStep._configure(self)
    self.dt = self.inventory.dt
    return


# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStepUniform.
  """
  return TimeStepUniform()


# End of file 
