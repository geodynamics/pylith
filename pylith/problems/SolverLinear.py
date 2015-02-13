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
    ModuleSolverLinear.__init__(self)
    return


  def initialize(self, fields, jacobian, formulation):
    """
    Initialize linear solver.
    """
    ModuleSolverLinear.initialize(self, fields, jacobian, formulation)
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
  return SolverLinear()


# End of file 
