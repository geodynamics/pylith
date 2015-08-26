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
## @brief Python PyLith abstract base class for progress monitor.
##
## Factory: progress_monitor

from pylith.utils.PetscComponent import PetscComponent
import datetime

# ProgressMonitor class
class ProgressMonitor(PetscComponent):
  """
  Python abstract base class for progress monitor.

  Factory: progress_monitor.
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
    ## @li \b update_percent Frequency of progress updates (percent).
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    updatePercent = pyre.inventory.float("update_percent", default=5.0)
    updatePercent.meta['tip'] = "Frequency of progress updates (percent)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="progressmonitor"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="progress_monitor")
    self.isMaster = True
    return


  def open(self):
    self.iupdate = None
    self.datetimeStart = datetime.datetime.now()
    
    try:
      import pylith.mpi.mpi as mpi
      self.isMaster = 0 == mpi.rank()
    except:
      self.isMaster = True
    if self.isMaster:
      self._open()
    return


  def close(self):
    if self.isMaster:
      self._close()
    return


  def update(self, current, start, stop):
    if not self.iupdate is None:
      percentComplete = (100*(current-start))/(stop-start)
    else:
      self.iupdate = 0
      percentComplete = 0.0
    if percentComplete >= self.iupdate*self.updatePercent:
      now = datetime.datetime.now()
      if percentComplete > 0.0:
        finished = self.datetimeStart + datetime.timedelta(seconds=100.0/percentComplete * ((now-self.datetimeStart).total_seconds()))
      else:
        finished = "TBD"
      if self.isMaster:
        self._update(current, start, stop, now, finished, percentComplete)
      self.iupdate += 1
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)

    self.updatePercent = self.inventory.updatePercent
    return

  def _open(self):
    return


  def _close(self):
    return


  def _update(self, current, start, stop, now, finished, percentComplete):
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
  """
  Factory associated with ProgressMonitor.
  """
  return ProgressMonitor()


# End of file 
