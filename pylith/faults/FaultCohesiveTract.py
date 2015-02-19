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

## @file pylith/faults/FaultCohesiveTract.py
##

## @brief Python object for a fault surface with dynamic
## (friction) fault implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveTract as ModuleFaultCohesiveTract

from pylith.utils.NullComponent import NullComponent

# FaultCohesiveTract class
class FaultCohesiveTract(FaultCohesive, Integrator, ModuleFaultCohesiveTract):
  """
  Python object for a fault surface with dynamic (friction) fault
  implemented with cohesive elements.

  Factory: fault
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivedyn"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoTr "

    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': ["normal_dir"],
            'data': ["slip",
                     "traction"]},
}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._info.log("Pre-initializing fault '%s'." % self.label())
    FaultCohesive.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)

    ModuleFaultCohesiveTract.quadrature(self, self.faultQuadrature)

    if mesh.dimension() == 2:
      self.availableFields['cell']['info'] += ["strike_dir"]
    elif mesh.dimension() == 3:
      self.availableFields['cell']['info'] += ["strike_dir",
                                               "dip_dir"]

    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    ModuleFaultCohesiveTract.verifyConfiguration(self, self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    self._info.log("Initializing fault '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    
    FaultCohesive.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    Integrator.poststep(self, t, dt, fields)
    FaultCohesive.poststep(self, t, dt, fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = ModuleFaultCohesiveTract.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveTract.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveTract.cellField(self, name)
    else:
      field = ModuleFaultCohesiveTract.cellField(self, name, fields)
    return field


  def finalize(self):
    """
    Cleanup.
    """
    FaultCohesive.finalize(self)
    Integrator.finalize(self)
    self.output.close()
    self.output.finalize()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveTract.
    """
    ModuleFaultCohesiveTract.__init__(self)
    return
    
  
  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFault("Fault", self)
    #self.perfLogger.logFields("Fault", self.fields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveTract.
  """
  return FaultCohesiveTract()


# End of file 
