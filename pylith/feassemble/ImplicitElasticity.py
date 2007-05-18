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

## @file pylith/feassemble/ImplicitElasticity.py
##
## @brief Python object for implicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorImplicit import IntegratorImplicit

# ImplicitElasticity class
class ImplicitElasticity(IntegratorImplicit):
  """
  Python object for implicit time integration of dynamic elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="implicitelasticity"):
    """
    Constructor.
    """
    IntegratorImplicit.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.ImplicitElasticity()
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
  
  
# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ImplicitElasticity.
  """
  return ImplicitElasticity()


# End of file 
