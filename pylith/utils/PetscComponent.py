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
  

  def compilePerformanceLog(self):
    """
    Compile performance and memory information.
    """
    for component in self.components():
      if isinstance(component, PetscComponent):
        component.compilePerformanceLog()
        if hasattr(component, 'perfLogger'):
          self.perfLogger.join(component.perfLogger)

      # Facility arrays are not PetscComponents but have components().
      elif hasattr(component, "components"):
        for subcomponent in component.components():
          if isinstance(subcomponent, PetscComponent):
            subcomponent.compilePerformanceLog()
            if hasattr(subcomponent, 'perfLogger'):
              self.perfLogger.join(subcomponent.perfLogger)
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
        for subcomponent in self.components():
          if isinstance(subcomponent, PetscComponent):
            subcomponent.cleanup()

    self._cleanup()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    return
    

# End of file 
