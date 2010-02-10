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

## @file pylith/problems/ExplicitLumpedTet4.py
##
## @brief Python ExplicitLumpedTet4 object for solving equations using an
## explicit formulation with a lumped Jacobian matrix that is stored
## as a Field.
##
## Factory: pde_formulation

from ExplicitLumped import ExplicitLumped
from pylith.utils.profiling import resourceUsageString

# ExplicitLumpedTet4 class
class ExplicitLumpedTet4(ExplicitLumped):
  """
  Python ExplicitLumpedTet4 object for solving equations using an explicit
  formulation.

  The formulation has the general form, [A(t)] {u(t+dt)} = {b(t)},
  where we want to solve for {u(t+dt)}, A(t) is usually constant
  (i.e., independent of time), and {b(t)} usually depends on {u(t)}
  and {u(t-dt)}.

  Jacobian: A(t)
  solution: u(t+dt)
  residual: b(t) - A(t) \hat u(t+dt)
  constant: b(t)

  Factory: pde_formulation.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicitlumpedtet4"):
    """
    Constructor.
    """
    ExplicitLumped.__init__(self, name)
    self._loggingPrefix = "TSEx "
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicitTet4 import ElasticityExplicitTet4
    return ElasticityExplicitTet4()


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    ExplicitLumped._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with ExplicitLumped.
  """
  return ExplicitLumpedTet4()


# End of file 
