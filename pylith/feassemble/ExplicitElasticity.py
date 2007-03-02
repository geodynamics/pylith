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

## @file pylith/feassemble/ExplicitElasticity.py
##
## @brief Python object for explicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorExplicit import IntegratorExplicit

# ExplicitElasticity class
class ExplicitElasticity(IntegratorExplicit):
  """
  Python object for explicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicitelasticity"):
    """
    Constructor.
    """
    IntegratorExplicit.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.ExplicitElasticity()
    return


  def initMaterial(self, mesh, material):
    """
    Initialize material properties.
    """
    self._info.log("Initializing integrator for material '%s'." % \
                   material.matname)
    material.initialize()
    
    self.material = material
    self.cppHandle.material = self.material.cppHandle
    #self.cppHandle.createParameters(mesh.cppHandle)
    return
  
  
# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ExplicitElasticity.
  """
  return ExplicitElasticity()


# End of file 
