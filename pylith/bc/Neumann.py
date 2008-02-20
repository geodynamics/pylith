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

## @file pylith/bc/Neumann.py
##
## @brief Python object for managing traction boundary conditions.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Integrator import Integrator

# Neumann class
class Neumann(BoundaryCondition, Integrator):
  """
  Python object for managing traction boundary conditions.

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
    ## @li None
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for numerical integration
    ## @li \b output Output manager associated with diagnostic output.

    import pyre.inventory

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."

    from pylith.meshio.OutputNeumann import OutputNeumann
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputNeumann)
    output.meta['tip'] = "Output manager associated with diagnostic output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="neumann"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "NeBC "
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["tractions"],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    self.quadrature.preinitialize()
    self.output.preinitialize(self)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    if self.quadrature.cellDim != self.mesh.dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.quadrature.cellDim,
               self.label, self.mesh.dimension()-1)    
    self.output.verifyConfiguration()

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize Neumann boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    self.cppHandle.quadrature = self.quadrature.cppHandle
    BoundaryCondition.initialize(self, totalTime, numTimeSteps)

    from pylith.topology.Mesh import Mesh
    self.boundaryMesh = Mesh()
    self.boundaryMesh.initialize(self.mesh.coordsys)
    self.cppHandle.boundaryMesh(self.boundaryMesh.cppHandle)

    if None != self.output:
      self.output.initialize(self.quadrature)
      self.output.writeInfo()

    self._logger.eventEnd(logEvent)
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.boundaryMesh, None, None)


  def getCellField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle)
    else:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle,
                                                    fields.cppHandle)
    return (field, fieldType)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.quadrature = self.inventory.quadrature
    self.output = self.inventory.output
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.bc.bc as bindings
      self.cppHandle = bindings.Neumann()    
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with Neumann.
  """
  return Neumann()

  
# End of file 
