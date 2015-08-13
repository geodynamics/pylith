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

## @file pylith/solver/ProgressMonitorTime.py
##
## @brief Python PyLith object for monitoring progress of time dependent problem.
##
## Factory: progress_monitor

from ProgressMonitor import ProgressMonitor

# ProgressMonitorTime class
class ProgressMonitorTime(ProgressMonitor):
  """
  Python PyLith object for monitoring progress of time dependent problem.

  Factory: progress_monitor.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ProgressMonitor.Inventory):
    """
    Python object for managing ProgressMonitorTime facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ProgressMonitorTime facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of output file.
    ## @li \b t_units Units for simulation time in output.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="progress.txt")
    filename.meta['tip'] = "Name of output file."

    tUnits = pyre.inventory.str("t_units", default="year")
    tUnits.meta['tip'] = "Units for simulation time in output."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="progressmonitortime"):
    """
    Constructor.
    """
    ProgressMonitor.__init__(self, name)
    self.fout = None
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    ProgressMonitor._configure(self)

    self.filename = self.inventory.filename
    self.tUnits = self.inventory.tUnits
    return


  def _open(self):
    import pyre.units
    uparser = pyre.units.parser()
    self.tSimScale = uparser.parse(self.tUnits)

    self.fout = open(self.filename, "w")
    self.fout.write("Timestamp                     Simulation t   % complete   Est. completion\n")
    self.fout.flush()
    return


  def _close(self):
    if self.fout:
      self.fout.close()
      self.fout = None
    return


  def _update(self, t, tStart, tEnd, now, finished, percentComplete):
    tSimNorm = t.value / self.tSimScale.value
    self.fout.write("%s   %8.2f*%s   %10.0f   %s\n" % (now, tSimNorm, self.tUnits, percentComplete, finished))
    self.fout.flush()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
  """
  Factory associated with ProgressMonitorTime.
  """
  return ProgressMonitorTime()


# End of file 
