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

## @file pyre/feassemble/Integrator.py
## @brief Python abstract base class for finite-element integration.

from pyre.components.Component import Component

# Integrator class
class Integrator(Component):
  """Python abstract base class for finite-element integration."""

  # PUBLIC METHODS /////////////////////////////////////////////////////

  # INITIALIZE
  # set mesh, basis, basisDer, quadrature

  def integrateResidual(self, element, state):
    """Integrate residual for element."""
    raise NotImplementedError, \
          "Integrator::integrateResidual() not implemented."
    return


  def integrateJacobian(self, element, state):
    """Integrate the Jacobian weak form over an element using the given
    quadrature."""
    raise NotImplementedError, \
          "Integrator::integrateJacobian() not implemented."
    return


  def __init__(self, name="integrator"):
    """Constructor."""
    Component.__init__(self, name, facility="integrator")
    self.mesh = None
    self.basis = None
    self.basisDer = None
    self.quadrature = None
    return


# version
__id__ = "$Id$"

# End of file 
