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
    return


  def initialize(self, mesh, field):
    """
    Initialize solver.
    """
    self._createCppHandle()
    self.cppHandle.initialize(mesh.cppHandle, field)
    return


  def solve(self, fieldOut, jacobian, fieldIn):
    """
    Solve linear system.
    """
    self._info.log("Solving linear equations.")
    assert(None != self.cppHandle)
    self.cppHandle.solve(fieldOut, jacobian, fieldIn)
    return
  

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Solver._configure(self)
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.solver.solver as bindings
      self.cppHandle = bindings.SolverLinear()
    return
  

# End of file 
