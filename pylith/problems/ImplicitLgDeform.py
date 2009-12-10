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
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with ImplicitLgDeform.
  """
  return ImplicitLgDeform()


# End of file 
