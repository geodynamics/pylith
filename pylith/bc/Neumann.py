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
from bc import Neumann as ModuleNeumann

# Neumann class
class Neumann(BoundaryCondition, Integrator, ModuleNeumann):
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

    from pylith.feassemble.Quadrature import SubMeshQuadrature
    quadrature = pyre.inventory.facility("quadrature",
                                         factory=SubMeshQuadrature)
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
    self.bcQuadrature.preinitialize(mesh.coordsys().spaceDim())
    self.quadrature(self.bcQuadrature)
    self.output.preinitialize(self)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self, self.mesh)
    Integrator.verifyConfiguration(self)
    if self.bcQuadrature.cellDim() != self.mesh.dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.bcQuadrature.cellDim(),
               self.label(), self.mesh.dimension()-1)    
    self.output.verifyConfiguration(self.mesh)
    ModuleNeumann.verifyConfiguration(self, self.mesh)

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize Neumann boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)
    
    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    #from pylith.topology.Mesh import Mesh
    #self.boundaryMesh = Mesh()
    #self.boundaryMesh.initialize(self.mesh.coordsys)
    #self.cppHandle.boundaryMesh(self.boundaryMesh.cppHandle)

    #if None != self.output:
    #  self.output.initialize(normalizer, self.quadrature)
    #  self.output.writeInfo()

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
      (field, fieldType) = self.cellField(name)
    else:
      (field, fieldType) = self.cellField(name, fields)
    return (field, fieldType)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.bcQuadrature = self.inventory.quadrature
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleNeumann.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with Neumann.
  """
  return Neumann()

  
# End of file 
