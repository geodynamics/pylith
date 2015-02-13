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

## @file pylith/problems/ExplicitLgDeform.py
##
## @brief Python ExplicitLgDeform object for solving equations using an
## explicit formulation with rigid body motion and small strains.
##
## Factory: pde_formulation

from Explicit import Explicit

# ExplicitLgDeform class
class ExplicitLgDeform(Explicit):
  """
  Python ExplicitLgDeform object for solving equations using an explicit
  formulation with rigid body motion and small strains.

  Factory: pde_formulation.
  """

  class Inventory(Explicit.Inventory):
    """
    Python object for managing ImplicitLgDeform facilities and properties.

    Provide appropriate solver for small strains as the default.
    """

    ## @class Inventory
    ## Python object for managing Explicit facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicitlgdeform"):
    """
    Constructor.
    """
    Explicit.__init__(self, name)
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicitLgDeform import ElasticityExplicitLgDeform
    return ElasticityExplicitLgDeform()


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Explicit._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with ExplicitLgDeform.
  """
  return ExplicitLgDeform()


# End of file 
