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

## @file pyre/feassemble/IntegratorElasticity.py
## @brief Python finite-element integrator for elasticity.

from Integrator import Integrator

# IntegratorElasticity class
class IntegratorElasticity(Integrator):
  """Python finite-element integrator for elasticity."""

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def integrateResidual(self, element, state):
    """Integrate residual for element."""

    (v0, jacobian, inv, detJ) = self.mesh.elementGeometry(element)

    return residual


  def integrateJacobian(self, element, state):
    """Integrate the Jacobian weak form over an element using the given
    quadrature."""

    (v0, jacobian, inv, detJ) = self.mesh.elementGeometry(element)

    # This is the inhomogeneous, anisotropic
    for q, point, weight in zip(range(len(self.quadrature.weights)),
                                self.quadrature.points,
                                self.quadrature.weights):
      u_q = 0.0
      for j in len(self.elementJac):
        u_q += u[j]*self.basis[q,j]
      D_q = material.elasticityConsts(point, u_q)
      for j in len(self.elementJac):
        t_der[0] = invJ[0]*self.basisDer[q,j,0] + invJ[2]*self.basisDer[q,j,1]
        t_der[1] = invJ[1]*self.basisDer[q,j,0] + invJ[3]*self.basisDer[q,j,1]
        for k in len(self.elementJac[0]):
          b_der[0] = invJ[0]*self.basisDer[q,k,0] + \
                     invJ[2]*self.basisDer[q,k,1]
          b_der[1] = invJ[1]*self.basisDer[q,k,0] + \
                     invJ[3]*self.basisDer[q,k,1]
          self.elementJac[j,k] += (t_der[0]*b_der[0] +
                                   t_der[1]*b_der[1])*D_q*weight*detJ
    return elementJac


  def __init__(self, name="integratorelasticity"):
    """Constructor."""
    Integrator.__init__(self, name)
    return


# version
__id__ = "$Id$"

# End of file 
