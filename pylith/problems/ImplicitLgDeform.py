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

## @file pylith/problems/ImplicitLgDeform.py
##
## @brief Python ImplicitLgDeform object for solving equations using
## an implicit formulation with rigid body motions and small strains.
##
## Factory: pde_formulation

from Implicit import Implicit

# ImplicitLgDeform class
class ImplicitLgDeform(Implicit):
  """
  Python ImplicitLgDeform object for solving equations using an implicit
  formulation with rigid body motions and small strains.

  Factory: pde_formulation.
  """

  class Inventory(Implicit.Inventory):
    """
    Python object for managing ImplicitLgDeform facilities and properties.

    Provide appropriate solver for small strains as the default.
    """

    ## @class Inventory
    ## Python object for managing ExplicitLumped facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    from SolverNonlinear import SolverNonlinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverNonlinear)
    solver.meta['tip'] = "Algebraic solver."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="implicitlgdeform"):
    """
    Constructor.
    """
    Implicit.__init__(self, name)
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityImplicitLgDeform import ElasticityImplicitLgDeform
    return ElasticityImplicitLgDeform()


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Implicit._configure(self)
    self.solver = self.inventory.solver
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with ImplicitLgDeform.
  """
  return ImplicitLgDeform()


# End of file 
