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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/ExplicitLumpedTri3.py
##
## @brief Python ExplicitLumpedTri3 object for solving equations using an
## explicit formulation with a lumped Jacobian matrix that is stored
## as a Field.
##
## Factory: pde_formulation

from ExplicitLumped import ExplicitLumped
from pylith.utils.profiling import resourceUsageString

# ExplicitLumpedTri3 class
class ExplicitLumpedTri3(ExplicitLumped):
  """
  Python ExplicitLumpedTri3 object for solving equations using an explicit
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

  def __init__(self, name="explicitlumpedtri3"):
    """
    Constructor.
    """
    ExplicitLumped.__init__(self, name)
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicitTri3 import ElasticityExplicitTri3
    integrator = ElasticityExplicitTri3()
    integrator.normViscosity(self.normViscosity)
    return integrator


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
  return ExplicitLumpedTri3()


# End of file 
