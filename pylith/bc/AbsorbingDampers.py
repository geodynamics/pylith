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

## @file pylith/bc/AbsorbingDampers.py
##
## @brief Python object for managing absorbing boundary condition
## using simple dashpots.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Integrator import Integrator
from bc import AbsorbingDampers as ModuleAbsorbingDampers

# AbsorbingDampers class
class AbsorbingDampers(BoundaryCondition, Integrator, ModuleAbsorbingDampers):
  """
  Python object for managing absorbing boundary condition using simple
  dashpots.

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

    import pyre.inventory

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="absorbingdampers"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    ModuleAbssorbingDampers.__init__(self)
    self._loggingPrefix = "AbBC "
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    self.bcQuadrature.preinitialize(mesh.coordsys().spaceDim())
    self.quadrature(self.bcQuadrature)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    if self.bcQuadrature.cellDim() != self.mesh.dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.bcQuadrature.cellDim,
               self.label, self.mesh.dimension()-1)    

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize AbsorbingDampers boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)    
    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    self._logger.eventEnd(logEvent)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.bcQuadrature = self.inventory.quadrature
    return


# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with AbsorbingDampers.
  """
  return AbsorbingDampers()

  
# End of file 
