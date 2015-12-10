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

## @file pylith/solver/Solver.py
##
## @brief Python PyLith abstract base class for solver.
##
## Factory: solver

from pylith.utils.PetscComponent import PetscComponent

# VALIDATORS ///////////////////////////////////////////////////////////

# Validate use of CUDA.
def validateUseCUDA(value):
  from pylith.utils.utils import isCUDAEnabled
  if value and not isCUDAEnabled:
    raise ValueError("PyLith is not built with CUDA support.")
  return value


# Solver class
class Solver(PetscComponent):
  """
  Python abstract base class for solver.

  Factory: solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Solver facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Solver facilities and properties.
    ##
    ## \b Properties
    ## @li \b use_cuda Use CUDA in solve if supported by solver.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    useCUDA = pyre.inventory.bool("use_cuda", default=False,
                                  validator=validateUseCUDA)
    useCUDA.meta['tip'] = "Enable use of CUDA for finite-element integrations."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solver"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solver")
    return


  def preinitialize(self):
    if self.useCUDA:
      # Set vec_type for CUDA, if it has not already been set.
      from pylith.utils.petsc import optionsSetValue, optionsHasName
      if not optionsHasName("-vec_type", 0):
        optionsSetValue("-vec_type", "cusp")
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)

    self.useCUDA = self.inventory.useCUDA
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return Solver()


# End of file 
