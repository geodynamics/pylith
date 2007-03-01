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

## @file pylith/utils/PetscManager.py

## @brief Python PetscManager object for managing PETSc options.

## The PetscManager also takes care of initializing and finalizing
## PETSc.

from pyre.components.Component import Component

# PetscManager class
class PetscManager(Component):
  """
  Python PetscManager object for managing PETSc options.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="petsc"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="petsc")
    self.options = []
    return


  def initialize(self):
    """
    Initialize PETSc.
    """
    self._info.log("Initializing PETSc.")
    import pylith.utils.petsc as bindings
    import sys
    args = [sys.executable]
    options = self._getOptions()
    if len(options) > 0:
      for arg in options:
        args.append(arg)
    bindings.petsc_initialize(args)
    return


  def finalize(self):
    """
    Finalize PETSc.
    """
    self._info.log("Finalizing PETSc.")
    import pylith.utils.petsc as bindings
    bindings.petsc_finalize()
    return
  

  def updateConfiguration(self, registry):
    """
    Update Pyre configuration.
    """
    self.options = [
      (name, descriptor.value) for name, descriptor in registry.properties.iteritems()
      ]
    return []


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _getOptions(self):
    """
    Cleanup options for PETSc.
    """
    args = []
    for iname, value in self.options:
      args.append('-' + iname)
      if value != 'true':
        args.append(value)
    return args


# End of file 
