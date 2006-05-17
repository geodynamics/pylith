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

## @file pyre/feassemble/Assembler.py
## @brief Python finite-element assembler.

from pyre.components.Component import Component

# Assembler class
class Assembler(Component):
  """Python finite-element assembler."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Assembler facilities and properties."""

    ## @class Inventory
    ## Python object for managing Assembler facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="assembler"):
    """Constructor."""
    Component.__init__(self, name, facility="assembler")
    return

  def integrateResidual(self, e):
    # This is v0, Jac, invJac, detJ
    geometry = self.mesh.elementGeometry(e)
    return self.integrateResidual(e, geometry)

  def integrateJacobian(self, u, e):
    # This is v0, Jac, invJac, detJ
    geometry = self.mesh.elementGeometry(e)
    return self.integrateJacobian(u, e, geometry)

  def integrateJacobian(self, u, e, geometry):
    '''Integrate the Jacobian weak form over an element e using the given quadrature'''
    det = geometry[3]
    Jinv = geometry[2]
    # This is the inhomogeneous, anisotropic
    for q, point, weight in zip(range(len(self.quadrature.weights)), self.quadrature.points, self.quadrature.weights):
      u_q = 0.0
      for j in len(self.elementJac):
        u_q += u[j]*self.basis[q,j]
      D_q = material.elasticityConstant(point, u_q)
      for j in len(self.elementJac):
        t_der[0] = Jinv[0]*self.basisDer[q,j,0] + Jinv[2]*self.basisDer[q,j,1]
        t_der[1] = Jinv[1]*self.basisDer[q,j,0] + Jinv[3]*self.basisDer[q,j,1]
        for k in len(self.elementJac[0]):
          b_der[0] = Jinv[0]*self.basisDer[q,k,0] + Jinv[2]*self.basisDer[q,k,1]
          b_der[1] = Jinv[1]*self.basisDer[q,k,0] + Jinv[3]*self.basisDer[q,k,1]
          self.elementJac[j,k] += (t_der[0]*b_der[0] + t_der[1]*b_der[1])*D_q*weight*det
    return self.elementJac

  def assembleResidual(self, res):
    for e in mesh.elements():
      res.update(e, self.integrateResidual(e), ADD_VALUES)
    return

  def assembleJacobian(self, jac, u):
    for e in mesh.elements():
      jac.updateOperator(field, e, self.integrateJacobian(u.restrict(e), e), ADD_VALUES)
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return
  

# version
__id__ = "$Id$"

# End of file 
