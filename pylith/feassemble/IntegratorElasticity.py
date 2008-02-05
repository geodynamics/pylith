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

    self._logger.eventEnd(logEvent)
    return
  
  
  def poststep(self, t, dt, totalTime):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.material.poststep(t, dt, totalTime)

    self._logger.eventEnd(logEvent)
    return


# End of file 
