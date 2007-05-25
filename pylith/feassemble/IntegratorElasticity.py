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

## @file pylith/feassemble/IntegratorElasticity.py
##
## @brief Python object implementing sgeneral methods for time
## integration of the elasticity equation using finite-elements.
##
## Factory: integrator

from Integrator import Integrator

# IntegratorElasticity class
class IntegratorElasticity(Integrator):
  """
  Python object implementing sgeneral methods for time integration of
  the elasticity equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorelasticity"):
    """
    Constructor.
    """
    Integrator.__init__(self)
    import journal
    self._info = journal.info(name)
    return


  def initMaterial(self, mesh, material):
    """
    Initialize material properties.
    """
    self._info.log("Initializing integrator material '%s'." % material.label)
    material.initialize(mesh)
    self.material = material
    self.cppHandle.material = self.material.cppHandle
    return
  
  
# End of file 
