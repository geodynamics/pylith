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

## @file pylith/solver/ProgressMonitorStep.py
##
## @brief Python PyLith object for monitoring progress of problem with steps.
##
## Factory: progress_monitor

from ProgressMonitor import ProgressMonitor

# ProgressMonitorStep class
class ProgressMonitorStep(ProgressMonitor):
  """
  Python PyLith object for monitoring progress of problem wit steps.

  Factory: progress_monitor.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ProgressMonitor.Inventory):
    """
    Python object for managing ProgressMonitorStep facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ProgressMonitorStep facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of output file.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="progress.txt")
    filename.meta['tip'] = "Name of output file."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="progressmonitorstep"):
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
    return


  def _open(self):
    self.fout = open(self.filename, "w")
    self.fout.write("Timestamp                     Step        % complete   Est. completion\n")
    self.fout.flush()
    return


  def _close(self):
    if self.fout:
      self.fout.close()
      self.fout = None
    return


  def _update(self, stepCurrent, stepStart, stepEnd, now, finished, percentComplete):
    self.fout.write("%s   %10d   %10.0f   %s\n" % (now, stepCurrent, percentComplete, finished))
    self.fout.flush()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def progress_monitor():
  """
  Factory associated with ProgressMonitorStep.
  """
  return ProgressMonitorStep()


# End of file 
