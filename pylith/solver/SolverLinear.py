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

# SolverLinear class
class SolverLinear(Solver):
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
    import pylith.solver.solver as bindings
    self.cppHandle = bindings.SolverLinear()
    return


  def initialize(self, mesh, field):
    """
    Initialize solver.
    """
    self.cppHandle.initialize(mesh.cppHandle, field)
    return


  def solve(self, fieldOut, jacobian, fieldIn):
    """
    Solve linear system.
    """
    self.cppHandle.solve(fieldOut, jacobian, fieldIn)
    return
  

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Solver._configure(self)
    return


# End of file 
