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

## @file pylith/problems/TimeStepUser.py
##
## @brief Python class for marching forward in time with user
## specified time steps.  Time step sizes are repeated until total
## time is reached.
##
## Format of file with user specified time steps:
##
## Notes: Whitespace is ignored. Comment lines begin with '//'.
##
## units = VALUE
##
## dt0
## dt1
## dt2
##
## Factory: time_step

from TimeStep import TimeStep

from pylith.utils.profiling import resourceUsageString

# TimeStepUser class
class TimeStepUser(TimeStep):
  """
  Python abstract base class for marching format in time with user
  specified time steps. Time step sizes are repeated until total time
  is reached.


  Format of file with user specified time steps:

  Notes: Whitespace is ignored. Comment lines begin with '//'.

  units = VALUE

  dt0
  dt1
  dt2

  Factory: time_step.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(TimeStep.Inventory):
    """
    Python object for managing TimeStepUser facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing TimeStepUser facilities and properties.
    ##
    ## \b Properties
    ## @li \b total_time Time duration for simulation.
    ## @li \b filename Name of file with time step sizes.
    ## @li \b loop_steps Loop over steps if true, otherwise keep
    ##   using last time step size.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    filename = pyre.inventory.str("filename", default="timesteps.txt")
    filename.meta['tip'] = "Name of file with tme step sizes."

    loopSteps = pyre.inventory.bool("loop_steps", default=False)
    loopSteps.meta['tip'] = "Loop over steps if true, otherwise keep " \
                            "using last time step size."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timestepuser"):
    """
    Constructor.
    """
    TimeStep.__init__(self, name)
    self.steps = []
    self.index = 0
    return


  def initialize(self):
    """
    Initialize time step algorithm.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    TimeStep.initialize(self)

    self._readSteps()
    assert(len(self.steps) > 0)
    assert(self.index == 0)
    self.dt = self.steps[self.index]

    self._logger.eventEnd(logEvent)
    return


  def numTimeSteps(self):
    """
    Get number of total time steps (or best guess if adaptive).
    """
    from pyre.units.time import second
    t = 0.0*second
    nsteps = 0
    index = 0
    while t <= self.totalTime:
      t += self.steps[index]
      index += 1
      if index >= len(self.steps):
        if self.loopSteps:
          index = 0
        else:
          index -= 1
      nsteps += 1
    return nsteps


  def timeStep(self, integrators):
    """
    Get time step for advancing forward in time.
    """
    self.dt = self.steps[self.index]
    if self.index+1 < len(self.steps):
      self.index += 1
    elif self.loopSteps:
      self.index = 0
    return self.dt

  
  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    TimeStep._configure(self)
    self.totalTime = self.inventory.totalTime
    self.filename = self.inventory.filename
    self.loopSteps = self.inventory.loopSteps
    return


  def _readSteps(self):
    """
    Read time step sizes from file.
    """
    fin = open(self.filename, "r")
    lines = fin.readlines()
    fin.close()

    from pyre.units.unitparser import parser
    unitparser = parser()

    self.steps = []
    for line in lines:
      fields = line.split()
      if line[0:2] != "//" and len(fields) > 0:
        if line[0:5] != "units":
          try:
            value = float(fields[0])
          except ValueError:
            raise IOError("Unable to time step specification '%s'.\n"
                          "Expected floating point value." % line.rstrip())
          self.steps.append(value*units)
        else:
          if len(fields) != 3:
            raise IOError("Unable to parse units specification.\n"
                          "Expected 'units = VALUE', got '%s'." % line.rstrip())
          try:
            units = unitparser.parse(fields[2])
          except NameError:
            raise IOError("Unable to parse units specification.\n"
                          "Expected 'units = VALUE', got '%s'." % line.rstrip())            
    if len(self.steps) == 0:
      raise IOError("No time step size values found in time step file '%s'." % \
                      self.filename)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def time_step():
  """
  Factory associated with TimeStepUser.
  """
  return TimeStepUser()


# End of file 
