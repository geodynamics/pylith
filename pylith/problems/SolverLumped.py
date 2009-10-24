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

## @file pylith/solver/SolverLumped.py

## @brief Python PyLith simple solver for system with a lumped (i.e.,
## diagonal) Jacobian matrix.

from Solver import Solver
from problems import SolverLumped as ModuleSolverLumped

# SolverLumped class
class SolverLumped(Solver, ModuleSolverLumped):
  """
  Python PyLith linear algebraic solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Solver.Inventory):
    """
    Python object for managing SolverLumped facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolverLumped facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solverlinear"):
    """
    Constructor.
    """
    Solver.__init__(self, name)
    ModuleSolverLumped.__init__(self)
    return


  def initialize(self, fields, jacobian, formulation):
    """
    Initialize linear solver.
    """
    ModuleSolverLumped.initialize(self, fields, jacobian, formulation)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Solver._configure(self)
    return


  def _cleanup(self):
    """
    Deallocate PETSc and local data structures.
    """
    self.deallocate()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return SolverLumped()


# End of file 
