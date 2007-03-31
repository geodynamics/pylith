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

## @file pylith/problems/Explicit.py
##
## @brief Python Explicit object for solving equations using an
## explicit formulation.
##
## Factory: pde_formulation

from Formulation import Formulation

# Explicit class
class Explicit(Formulation):
  """
  Python Explicit object for solving equations using an explicit
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

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Explicit facilities and properties.
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

  def __init__(self, name="explicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    return


  def initialize(self, mesh, materialsBin, spaceDim):
    """
    Create explicit integrators for each element family.
    """
    from pylith.feassemble.ExplicitElasticity import ExplicitElasticity
    
    self._info.log("Initializing integrators.")
    self.integrators = []
    for material in materialsBin.materials:
      if material.quadrature.spaceDim != spaceDim:
        raise ValueError, \
              "Spatial dimension of problem is '%d' but quadrature " \
              "for material '%s' is for spatial dimension '%d'." % \
              (spaceDim, material.label, material.quadrature.spaceDim)
      integrator = ExplicitElasticity()
      integrator.setMesh(mesh)
      integrator.initQuadrature(material.quadrature)
      integrator.initMaterial(mesh, material)
      self.integrators.append(integrator)

    self._info.log("Creating fields and matrices.")
    self.dispT = mesh.cppHandle.createRealSection("dispT", spaceDim)
    self.dispTmdt = mesh.cppHandle.createRealSection("dispTmdt", spaceDim)
    self.dispTpdt = mesh.cppHandle.createRealSection("dispTpdt", spaceDim)
    self.constant = mesh.cppHandle.createRealSection("constant", spaceDim)
    self.jacobian = mesh.cppHandle.createMatrix(self.constant)

    self._info.log("Integrating Jacobian of operator.")
    for integrator in self.integrators:
      integrator.integrateJacobian(self.jacobian, self.dispT)
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Explicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: Explicit::prestep() not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating constant term in operator.")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(self.constant)
    for integrator in self.integrators:
      integrator.integrateConstant(self.constant, self.dispT, self.dispTmdt)

    self._info.log("Solving equations.")
    # solve
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    tmp = self.dispTmdt
    self.dispTmdt = self.dispT
    self.dispT = self.dispTpdt
    self.dispTpdt = tmp
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


  def _calcConstant(self):
    """
    Compute residual, {b(t)}.
    """
    self._info.log("WARNING: Explicit::calcConstant() not implemented.")
    return


  def _calcJacobian(self):
    """
    Compute Jacobian, [A(t)].
    """
    self._info.log("WARNING: Explicit::calcJacobian() not implemented.")
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
