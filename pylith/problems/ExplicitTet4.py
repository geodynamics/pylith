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

## @file pylith/problems/ExplicitTet4.py
##
## @brief Python ExplicitTet4 object for solving equations using an
## explicit formulation with a lumped Jacobian matrix that is stored
## as a Field.
##
## Factory: pde_formulation

from Explicit import Explicit
from pylith.utils.profiling import resourceUsageString

# ExplicitTet4 class
class ExplicitTet4(Explicit):
  """
  Python ExplicitTet4 object for solving equations using an explicit
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

  def __init__(self, name="explicittet4"):
    """
    Constructor.
    """
    Explicit.__init__(self, name)
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicitTet4 import ElasticityExplicitTet4
    integrator = ElasticityExplicitTet4()
    integrator.normViscosity(self.normViscosity)
    return integrator


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
  Factory associated with Explicit.
  """
  return ExplicitTet4()


# End of file 
