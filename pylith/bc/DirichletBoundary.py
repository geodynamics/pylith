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

from DirichletBC import DirichletBC
from DirichletBC import validateDOF
from bc import DirichletBoundary as ModuleDirichletBoundary

# DirichletBoundary class
class DirichletBoundary(DirichletBC, ModuleDirichletBoundary):
  """
  Python object for managing a DirichletBoundary (prescribed displacements)
  boundary condition with points on a surface.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(DirichletBC.Inventory):
    """
    Python object for managing DirichletBoundary facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing DirichletBoundary facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b output Output manager associated with diagnostic output.

    import pyre.inventory

    from pylith.meshio.OutputDirichlet import OutputDirichlet
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputDirichlet)
    output.meta['tip'] = "Output manager associated with diagnostic output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichletboundary"):
    """
    Constructor.
    """
    DirichletBC.__init__(self, name)
    self._loggingPrefix = "DiBC "
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
    DirichletBC.preinitialize(self, mesh)
    self.output.preinitialize(self)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    DirichletBC.verifyConfiguration(self)
    self.output.verifyConfiguration(self.mesh)

    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize DirichletBoundary boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    DirichletBC.initialize(self, totalTime, numTimeSteps, normalizer)

    self.output.initialize(normalizer)
    self.output.writeInfo()

    self._logger.eventEnd(logEvent)    
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    label = ""
    labelId = 0
    return (self.boundaryMesh(), label, labelId)


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = self.vertexField(name)
    else:
      field = self.vertexField(name, fields)
    return field


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    DirichletBC._configure(self)
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleDirichletBoundary.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with DirichletBoundary.
  """
  return DirichletBoundary()

  
# End of file 
