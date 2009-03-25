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

## @file pylith/solver/SolverLinear.py

## @brief Python PyLith linear algebraic solver.

from Solver import Solver
from problems import SolverLinear as ModuleSolverLinear

# SolverLinear class
class SolverLinear(Solver, ModuleSolverLinear):
  """
  Python PyLith linear algebraic solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Solver.Inventory):
    """
    Python object for managing SolverLinear facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolverLinear facilities and properties.
    ##
    ## \b Properties
    ## @li \b initial_guess_zero Use zero for initial guess.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    guessZero = pyre.inventory.bool("initial_guess_zero", default=True)
    guessZero.meta['tip'] = "Use zero for initial guess."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solverlinear"):
    """
    Constructor.
    """
    Solver.__init__(self, name)
    ModuleSolverLinear.__init__(self)
    return


  def initialize(self, fields, jacobian, formulation):
    """
    Initialize linear solver.
    """
    ModuleSolverLinear.initialize(self, fields, jacobian, formulation)
    self.initialGuessZero(self.guessZero)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Solver._configure(self)
    self.guessZero = self.inventory.guessZero
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return SolverLinear()


# End of file 
