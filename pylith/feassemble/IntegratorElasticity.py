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
    self.mesh = mesh
    self.output = material.output
    self.availableFields = material.availableFields
    self.materialObj = material

    Integrator.preinitialize(self, mesh)
    material.preinitialize(mesh)
    #self.output.preinitialize(self)

    # Set integrator's quadrature using quadrature from material
    self.quadrature(material.quadrature)
    self.material(material)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    Integrator.verifyConfiguration(self)
    #self.output.verifyConfiguration(self.mesh)

    self._logger.eventEnd(logEvent)    
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize material properties.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Initializing integrator for material '%s'." % \
                   self.materialObj.label)

    print "DDD"
    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    print "EEE"

    #self.output.initialize(normalizer, self.materialObj.quadrature)
    #self.output.writeInfo()
    #self.output.open(totalTime, numTimeSteps)

    self._logger.eventEnd(logEvent)
    return
  
  
  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    Integrator.poststep(self, t, dt, totalTime, fields)

    self._info.log("Writing material data.")
    #self.output.writeData(t+dt, fields)

    self._logger.eventEnd(logEvent)
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
      (field, fieldType) = self.cellField(name)
    else:
      (field, fieldType) = self.cellField(name, fields)
    return (field, fieldType)


# End of file 
