#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
    ## @li \b db Database of boundary condition parameters

    import pyre.inventory

    from pylith.feassemble.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", factory=SimpleDB, 
                                 family="spatial_database")
    db.meta['tip'] = "Database of boundary condition parameters."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="absorbingdampers"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
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
    self.createSubMesh(mesh)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self, self.mesh())
    Integrator.verifyConfiguration(self)
    if self.bcQuadrature.cellDim() != self.mesh().dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.bcQuadrature.cellDim,
               self.label(), self.mesh().dimension()-1)    
    ModuleAbsorbingDampers.verifyConfiguration(self, self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize AbsorbingDampers boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing absorbing boundary '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)    
    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return
  

  def finalize(self):
    """
    Cleanup.
    """
    BoundaryCondition.finalize(self)
    Integrator.finalize(self)
    self._modelMemoryUse()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    self.bcQuadrature = self.inventory.quadrature
    ModuleAbsorbingDampers.db(self, self.inventory.db)
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleAbsorbingDampers.__init__(self)
    return
  

  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFields("BoundaryConditions", self.parameterFields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with AbsorbingDampers.
  """
  return AbsorbingDampers()

  
# End of file 
