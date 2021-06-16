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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/faults/FaultCohesiveImpulses.py
##

## @brief Python object for a fault surface with slip impulses for
## Green's function implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveImpulses as ModuleFaultCohesiveImpulses

def validateDOF(value):
  """
  Validate list of fixed degrees of freedom.
  """
  try:
    size = len(value)
    num = map(int, value)
    for v in num:
      if v < 0:
        raise ValueError
  except:
    raise ValueError, \
          "'impuluse_dof' must be a zero based list of indices of degrees of " \
          "freedom at a vertex."
  return num
  

# FaultCohesiveImpulses class
class FaultCohesiveImpulses(FaultCohesive, Integrator, ModuleFaultCohesiveImpulses):
  """
  Python object for a fault surface with slip impulses for Green's
  functions implemented with cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesiveImpulses facilities and properties.
  
  \b Properties
  @li \b threshold Threshold for non-zero amplitude.
  @li \b impulse_dof Indices of slip components for impulses.
  
  \b Facilities
  @li \b db_impulse_amplitude Amplitude of slip impulses.
  @li \b output Output manager associated with fault data.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  from pyre.units.length import m

  threshold = pyre.inventory.dimensional("threshold", default=1.0e-6*m)
  threshold.meta['tip'] = "Threshold for non-zero amplitude."

  impulseDOF = pyre.inventory.list("impulse_dof", default=[], validator=validateDOF)
  impulseDOF.meta['tip'] = "Indices of impulse components " \
      "(0=1st DOF, 1=2nd DOF, etc)."

  from spatialdata.spatialdb.SimpleDB import SimpleDB
  dbImpulseAmp = pyre.inventory.facility("db_impulse_amplitude", family="spatial_database", factory=SimpleDB)
  dbImpulseAmp.meta['tip'] = "Amplitude of slip impulses."
  
  from pylith.meshio.OutputFaultImpulses import OutputFaultImpulses
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputFaultImpulses)
  output.meta['tip'] = "Output manager associated with fault data."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesiveimpulses"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoIm "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal_dir",
                     "impulse_amplitude",
                     "area",],
            'data': ["slip",
                     "traction_change"]},
         'cell': \
           {'info': ["partition"],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Pre-initializing fault '%s'." % self.label())
    FaultCohesive.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)

    ModuleFaultCohesiveImpulses.quadrature(self, self.faultQuadrature)

    if mesh.dimension() == 2:
      self.availableFields['vertex']['info'] += ["strike_dir"]
    elif mesh.dimension() == 3:
      self.availableFields['vertex']['info'] += ["strike_dir",
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
    ModuleFaultCohesiveImpulses.verifyConfiguration(self, self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize cohesive elements.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing fault '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    FaultCohesive.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = ModuleFaultCohesiveImpulses.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveImpulses.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveImpulses.cellField(self, name)
    else:
      field = ModuleFaultCohesiveImpulses.cellField(self, name, fields)
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
    import numpy
    FaultCohesive._configure(self)
    self.output = self.inventory.output

    ModuleFaultCohesiveImpulses.threshold(self, self.inventory.threshold.value)
    impulseDOF = numpy.array(self.inventory.impulseDOF, dtype=numpy.int32)
    ModuleFaultCohesiveImpulses.impulseDOF(self, impulseDOF)
    ModuleFaultCohesiveImpulses.dbImpulseAmp(self, self.inventory.dbImpulseAmp)
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveImpulses.
    """
    ModuleFaultCohesiveImpulses.__init__(self)
    return
    
  
  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFault("Fault", self)
    self.perfLogger.logFields("Fault", self.fields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveImpulses.
  """
  return FaultCohesiveImpulses()


# End of file 
