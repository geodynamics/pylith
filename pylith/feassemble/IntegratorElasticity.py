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

## @file pylith/feassemble/IntegratorElasticity.py
##
## @brief Python object implementing sgeneral methods for time
## integration of the elasticity equation using finite-elements.
##
## Factory: integrator

from Integrator import Integrator

# IntegratorElasticity class
class IntegratorElasticity(Integrator):
  """
  Python object implementing sgeneral methods for time integration of
  the elasticity equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorelasticity"):
    """
    Constructor.
    """
    Integrator.__init__(self)
    self.output = None
    self.availableFields = None
    self.name = "Integrator Elasticity"

    # Setup journal (not a Component, so not setup already)
    import journal
    self._info = journal.info(name)
    return


  def preinitialize(self, mesh, material):
    """
    Setup integrator.
    """
    import weakref
    self.mesh = weakref.ref(mesh)
    self.output = material.output
    self.availableFields = material.availableFields
    self.materialObj = material

    Integrator.preinitialize(self, mesh)
    material.preinitialize(mesh)
    self.output.preinitialize(self)

    # Set integrator's quadrature using quadrature from material
    self.quadrature(material.quadrature)
    self.material(material)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    Integrator.verifyConfiguration(self)
    self.materialObj.verifyConfiguration()

    if self.mesh().dimension() != self.materialObj.dimension():
      raise ValueError("Mesh dimension is '%d' but material '%s' of type " \
                         "'%s' applies to dimension '%d'." % \
                       (self.mesh().dimension(),
                        self.materialObj.label(),
                        self.materialObj,
                        self.materialObj.dimension()))
    self._verifyConfiguration()
    self.output.verifyConfiguration(self.mesh())

    self._eventLogger.eventEnd(logEvent)    
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize material properties.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    if 0 == comm.rank:
      self._info.log("Initializing integrator for material '%s'." % \
                       self.materialObj.label())
    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    self.materialObj.normalizer(normalizer)

    self._eventLogger.eventEnd(logEvent)
    return
  
  
  def writeData(self, t, fields):
    """
    Hook for writing data at time t.
    """
    logEvent = "%swrite" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    if 0 == comm.rank:
      self._info.log("Writing material data.")
    self.output.writeData(t, fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def finalize(self):
    """
    Cleanup.
    """
    Integrator.finalize(self)
    self.materialObj.finalize()
    self.output.close()
    self.output.finalize()
    self._modelMemoryUse()
    return


  def cleanup(self):
    """
    Deallocate locally managed data structures.
    """
    self.deallocate()
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return self.materialObj.getDataMesh()


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = self.cellField(name, self.mesh())
    else:
      field = self.cellField(name, self.mesh(), fields)
    return field


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initializeOutput(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize output.
    """
    self.output.initialize(normalizer, self.materialObj.quadrature)
    self.output.writeInfo()
    self.output.open(totalTime, numTimeSteps)
    return


  def _verifyConfiguration(self):
    raise NotImplementedError("Implement _verifyConfiguration() in child class.")


  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.materialObj.perfLogger.logFields("Output", self.outputFields())
    return


# End of file 
