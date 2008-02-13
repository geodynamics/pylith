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
    import journal
    self._info = journal.info(name)
    self.output = None
    self.availableFields = None
    self.name = "Integrator Elasticity"
    return


  def preinitialize(self, mesh, material):
    """
    Setup integrator.
    """
    Integrator.preinitialize(self, mesh)
    
    assert(None != self.cppHandle)
    self.mesh = mesh

    material.preinitialize()

    self.quadrature = material.quadrature
    self.cppHandle.quadrature = self.quadrature.cppHandle

    self.material = material
    self.cppHandle.material = self.material.cppHandle
    self.output = material.output
    self.availableFields = material.availableFields
    self.output.preinitialize(self)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    Integrator.verifyConfiguration(self)
    self.material.verifyConfiguration()
    self.output.verifyConfiguration()    
    return


  def initialize(self, totalTime, numTimeSteps):
    """
    Initialize material properties.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Initializing integrator for material '%s'." % \
                   self.material.label)

    self.material.initialize(self.mesh, totalTime, numTimeSteps)
    self.output.initialize(self.quadrature.cppHandle)
    self.output.writeInfo()
    self.output.open(totalTime, numTimeSteps)

    self._logger.eventEnd(logEvent)
    return
  
  
  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Writing material data.")
    self.output.writeData(t+dt, fields)

    self._logger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return self.material.getDataMesh()


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle)
    else:
      (field, fieldType) = self.cppHandle.cellField(name, self.mesh.cppHandle,
                                                   fields.cppHandle)
    return (field, fieldType)


# End of file 
