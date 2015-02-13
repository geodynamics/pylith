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

## @file pylith/solver/SolverNonlinear.py

## @brief Python PyLith nonlinear algebraic solver.

from Solver import Solver
from problems import SolverNonlinear as ModuleSolverNonlinear

# SolverNonlinear class
class SolverNonlinear(Solver, ModuleSolverNonlinear):
  """
  Python PyLith nonlinear algebraic solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Solver.Inventory):
    """
    Python object for managing SolverNonlinear facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolverNonlinear facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solvernonlinear"):
    """
    Constructor.
    """
    Solver.__init__(self, name)
    ModuleSolverNonlinear.__init__(self)
    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Solver._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return SolverNonlinear()


# End of file 
