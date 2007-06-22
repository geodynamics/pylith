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

## @file pylith/problems/Formulation.py
##
## @brief Python abstract base class for formulations of solving equations.
##
## Factory: pde_formulation

from pyre.components.Component import Component

# Formulation class
class Formulation(Component):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Formulation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.
    ## @li \b output Solution output.

    import pyre.inventory

    from pylith.solver.SolverLinear import SolverLinear
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLinear)
    solver.meta['tip'] = "Algebraic solver."

    from pylith.meshio.SingleOutput import SingleOutput
    output = pyre.inventory.facility("output", family="object_bin",
                                     factory=SingleOutput)
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="pde_formulation")
    self.integrators = None
    self.constraints = None
    self.fields = None
    self.solnField = None
    self._istep = 0
    return


  def initialize(self, mesh, materials, boundaryConditions,
                 interfaceConditions, dimension, dt):
    """
    Create integrators for each element family.
    """
    from pylith.feassemble.Integrator import implementsIntegrator
    from pylith.feassemble.Constraint import implementsConstraint

    from pylith.topology.FieldsManager import FieldsManager
    self.fields = FieldsManager(mesh)
    self.integrators = []
    self.constraints = []

    self._info.log("Initializing materials.")
    for material in materials.bin:
      if material.quadrature.spaceDim != dimension:
        raise ValueError, \
              "Spatial dimension of problem is '%d' but quadrature " \
              "for material '%s' is for spatial dimension '%d'." % \
              (dimension, material.label, material.quadrature.spaceDim)
      integrator = self.elasticityIntegrator()
      if not implementsIntegrator(integrator):
        raise TypeError, \
              "Could not use '%s' as an integrator for material '%s'. " \
              "Functionality missing." % (integrator.name, material.label)
      integrator.setMesh(mesh)
      integrator.initQuadrature(material.quadrature)
      integrator.initMaterial(mesh, material)
      self.integrators.append(integrator)
      self._info.log("Added elasticity integrator for material '%s'." % \
                     material.label)

    self._info.log("Initializing boundary conditions.")
    for bc in boundaryConditions.bin:
      bc.initialize(mesh)
      foundType = False
      if implementsIntegrator(bc):
        foundType = True
        self.integrators.append(bc)
        self._info.log("Added boundary condition '%s' as an integrator." % \
                       bc.label)
      if implementsConstraint(bc):
        foundType = True
        self.constraints.append(bc)
        self._info.log("Added boundary condition '%s' as a constraint." % \
                       bc.label)
      if not foundType:
        raise TypeError, \
              "Could not determine whether boundary condition '%s' is an " \
              "integrator or a constraint." % bc.name

    self._info.log("Initializing interior interfaces.")
    for ic in interfaceConditions.bin:
      ic.initialize(mesh)
      foundType = False
      if implementsIntegrator(ic):
        foundType = True
        self.integrators.append(ic)
        self._info.log("Added interface condition '%s' as an integrator." % \
                       ic.label)
      if implementsConstraint(ic):
        foundType = True
        self.constraints.append(ic)
        self._info.log("Added interface condition '%s' as a constraint." % \
                       ic.label)
      if not foundType:
        raise TypeError, \
              "Could not determine whether interface condition '%s' is an " \
              "integrator or a constraint." % ic.name

    self._info.log("Setting up solution output.")
    for output in self.output.bin:
      output.open(mesh)
      output.writeTopology()

    self._info.log("Creating solution field.")
    solnName = self.solnField['name']
    self.fields.addReal(solnName)
    self.fields.solutionField(solnName)
    self.fields.setFiberDimension(solnName, dimension)
    for constraint in self.constraints:
      constraint.setConstraintSizes(self.fields.getSolution())
    self.fields.allocate(solnName)
    for constraint in self.constraints:
      constraint.setConstraints(self.fields.getSolution())
    return


  def poststep(self, t, dt, totalTime):
    """
    Hook for doing stuff after advancing time step.
    """
    field = self.fields.getSolution()
    for output in self.output.bin:
      output.writeField(t+dt, self._istep, field, self.solnField['label'])
    self._istep += 1
    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    for integrator in self.integrators:
      integrator.finalize()
    for constraint in self.constraints:
      constraint.finalize()
    for output in self.output.bin:
      output.close()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.solver = self.inventory.solver
    self.output = self.inventory.output
    return


  def _reformJacobian(self, t, dt):
    """
    Reform Jacobian matrix for operator.
    """
    self._info.log("Reforming Jacobian of operator.")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, t+dt, self.fields)
    petsc.mat_assemble(self.jacobian)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Formulation.
  """
  return Formulation()


# End of file 
