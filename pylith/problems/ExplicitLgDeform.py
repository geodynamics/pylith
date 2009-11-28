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
