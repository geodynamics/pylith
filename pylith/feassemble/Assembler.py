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

  # PUBLIC METHODS /////////////////////////////////////////////////////

  # INITIALIZE
  # Set mesh, chart, integrator
  
  def assembleResidual(self, residual, state):
    """Assemble residual for field."""
    elements = self.mesh.restrict(self.chart)
    for element in elements:
      stateElem = self.mesh.restrict(element)
      residual.update(element,
                      self.integrator.integrateResidual(element, stateElem),
                      ADD_VALUES)
    return


  def assembleJacobian(self, jacobian, state):
    elements = self.mesh.restrict(self.chart)
    for element in elements:
      stateElem = self.mesh.restrict(element)
      jacobian.updateOperator(field, element,
                              self.integrator.integrateJacobian(element, stateElem),
                              ADD_VALUES)
    return


  def __init__(self, name="assembler"):
    """Constructor."""
    Component.__init__(self, name, facility="assembler")
    self.mesh = None
    self.chart = None
    self.integrator = None
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """Set members based using inventory."""
    return
  

# version
__id__ = "$Id$"

# End of file 
