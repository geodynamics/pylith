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

## @file pylith/utils/PetscComponent.py
##
## @brief Python PetscComponent object for aid in deallocating data
## structures before calling PetscFinalize().

from pyre.components.Component import Component

# PetscComponent class
class PetscComponent(Component):
  """
  Python PetscComponent object for aid in deallocating data structures
  before calling PetscFinalize().
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name, facility):
    """
    Constructor.
    """
    Component.__init__(self, name, facility)
    return
  

  def compilePerformanceLog(self, parentLogger):
    """
    Compile performance and memory information.
    """
    if hasattr(self, 'perfLogger'):
      if not parentLogger is None:
        parentLogger.join(self.perfLogger)

    for component in self.components():
      if isinstance(component, PetscComponent):
        component.compilePerformanceLog(parentLogger)

      # Facility arrays are not PetscComponents but have components().
      elif hasattr(component, "components"):
        for subcomponent in component.components():
          if isinstance(subcomponent, PetscComponent):
            subcomponent.compilePerformanceLog(parentLogger)
    return


  def cleanup(self):
    """
    Deallocate data structures.
    """
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.cleanup()

      # Facility arrays are not PetscComponents but have components().
      elif hasattr(component, "components"):
        for subcomponent in component.components():
          if isinstance(subcomponent, PetscComponent):
            subcomponent.cleanup()

    self._cleanup()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    if "deallocate" in dir(self):
      self.deallocate()
    return
    

# End of file 
