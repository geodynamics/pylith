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

## @file pylith/bc/DirichletBoundary.py
##
## @brief Python object for managing a Dirichlet (prescribed
## displacements) boundary condition with points on a surface.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Constraint import Constraint

def validateDOF(value):
  """
  Validate list of fixed degrees of freedom.
  """
  try:
    size = len(value)
    num = map(int, value)
    for v in num:
      if v < 0:
        raise ValueError
  except:
    raise ValueError, \
          "'fixed_dof' must be a zero based list of indices of fixed " \
          "degrees of freedom."
  return num
  

# DirichletBoundary class
class DirichletBoundary(BoundaryCondition, Constraint):
  """
  Python object for managing a DirichletBoundary (prescribed displacements)
  boundary condition with points on a surface.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryCondition.Inventory):
    """
    Python object for managing BoundaryCondition facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BoundaryCondition facilities and properties.
    ##
    ## \b Properties
    ## @li \b fixed_dof Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc).
    ## @li \b reference_t Reference time for rate of change of values.
    ##
    ## \b Facilities
    ## @li \b initial_db Database of parameters for initial values.
    ## @li \b rate_db Database of parameters for rate of change of values.
    ## @li \b output Output manager associated with diagnostic output.

    import pyre.inventory

    fixedDOF = pyre.inventory.list("fixed_dof", default=[],
                                   validator=validateDOF)
    fixedDOF.meta['tip'] = "Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc)."

    from pyre.units.time import s
    tRef = pyre.inventory.dimensional("reference_t", default=0.0*s)
    tRef.meta['tip'] = "Reference time for rate of change of values."

    from FixedDOFDB import FixedDOFDB
    db = pyre.inventory.facility("db", factory=FixedDOFDB,
                                 family="spatial_database")
    db.meta['tip'] = "Database of parameters for initial values."

    dbRate = pyre.inventory.facility("rate_db", factory=FixedDOFDB,
                                 family="spatial_database")
    dbRate.meta['tip'] = "Database of parameters for rate of change of values."

    from pylith.meshio.OutputDirichlet import OutputDirichlet
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputDirichlet)
    output.meta['tip'] = "Output manager associated with diagnostic output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichletpoints"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Constraint.__init__(self)
    self._loggingPrefix = "DiBC "
    self.fixedDOF = []
    self.availableFields = \
        {'vertex': \
           {'info': ["initial", "rate-of-change"],
            'data': []},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Constraint.preinitialize(self, mesh)
    self.cppHandle.fixedDOF = self.fixedDOF
    self.output.preinitialize(self)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self)
    Constraint.verifyConfiguration(self)
    self.output.verifyConfiguration(self.mesh)

    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize DirichletBoundary boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    assert(None != self.cppHandle)
    self.cppHandle.referenceTime = self.tRef.value
    self.dbRate.initialize()
    self.cppHandle.dbRate = self.dbRate.cppHandle

    BoundaryCondition.initialize(self, totalTime, numTimeSteps)

    from pylith.topology.Mesh import Mesh
    self.boundaryMesh = Mesh()
    self.boundaryMesh.initialize(self.mesh.coordsys)
    self.cppHandle.boundaryMesh(self.boundaryMesh.cppHandle)

    self.output.initialize()
    self.output.writeInfo()

    self._logger.eventEnd(logEvent)    
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.boundaryMesh, None, None)


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      (field, fieldType) = self.cppHandle.vertexField(name,
                                                      self.mesh.cppHandle)
    else:
      (field, fieldType) = self.cppHandle.vertexField(name,
                                                     self.mesh.cppHandle,
                                                     fields.cppHandle)
    return (field, fieldType)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.tRef = self.inventory.tRef
    self.fixedDOF = self.inventory.fixedDOF
    self.dbRate = self.inventory.dbRate
    self.output = self.inventory.output
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.bc.bc as bindings
      self.cppHandle = bindings.DirichletBoundary()    
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with DirichletBoundary.
  """
  return DirichletBoundary()

  
# End of file 
