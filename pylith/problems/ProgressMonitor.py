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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/solver/ProgressMonitor.py
##
## @brief Python PyLith abstract base class for solver.
##
## Factory: solver

from pylith.utils.PetscComponent import PetscComponent
import datetime

# ProgressMonitor class
class ProgressMonitor(PetscComponent):
  """
  Python abstract base class for solver.

  Factory: solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing ProgressMonitor facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ProgressMonitor facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of output file.
    ## @li \b t_units Units for simulation time in output.
    ## @li \b update_percent Frequency of progress updates (percent).
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="progress.txt")
    filename.meta['tip'] = "Name of output file."

    tUnits = pyre.inventory.str("t_units", default="year")
    tUnits.meta['tip'] = "Units for simulation time in output."

    updatePercent = pyre.inventory.float("update_percent", default=5.0)
    updatePercent.meta['tip'] = "Frequency of progress updates (percent)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="progress_monitor"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="progress_monitor")
    self.fout = None
    return


  def open(self):
    self.tPrev = None
    self.datetimeStart = datetime.datetime.now()
    
    import pyre.units
    uparser = pyre.units.parser()
    self.tSimScale = uparser.parse(self.tUnits)

    try:
      import pylith.mpi.mpi as mpi
      self.isMaster = 0 == mpi.rank()
    except:
      self.isMaster = True
    if self.isMaster:
      self.fout = open(self.filename, "w")
      self.fout.write("Timestamp                     Simulation t   % complete   Est. completion\n")
    return


  def close(self):
    if self.fout:
      self.fout.close()
      self.fout = None
    return


  def update(self, t, tStart, tEnd):
    writeUpdate = False
    if self.tPrev:
      incrCompleted = (t-self.tPrev) / (tEnd-tStart)
    else:
      incrCompleted = 0.0
    if not self.tPrev or incrCompleted > self.updatePercent/100.0:
      percentComplete = (t-tStart)/(tEnd-tStart)*100.0
      tSimNorm = t.value / self.tSimScale.value
      now = datetime.datetime.now()
      if percentComplete > 0.0:
        finished = self.datetimeStart + datetime.timedelta(seconds=100.0/percentComplete * ((now-self.datetimeStart).total_seconds()))
      else:
        finished = "TBD"
      if self.isMaster:
        self.fout.write("%s   %8.2f*%s   %10.0f   %s\n" % (now, tSimNorm, self.tUnits, percentComplete, finished))
      self.tPrev = t
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)

    self.filename = self.inventory.filename
    self.tUnits = self.inventory.tUnits
    self.updatePercent = self.inventory.updatePercent
    return


# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
  """
  Factory associated with ProgressMonitor.
  """
  return ProgressMonitor()


# End of file 
