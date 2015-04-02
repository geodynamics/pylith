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

## @file pylith/faults/FaultCohesiveKin.py
##

## @brief Python object for a fault surface with kinematic
## (prescribed) slip implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveKin as ModuleFaultCohesiveKin

# ITEM FACTORIES ///////////////////////////////////////////////////////

def eqsrcFactory(name):
  """
  Factory for earthquake source items.
  """
  from pyre.inventory import facility
  from EqKinSrc import EqKinSrc
  return facility(name, family="eq_kinematic_src", factory=EqKinSrc)


# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive, Integrator, ModuleFaultCohesiveKin):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesiveKin facilities and properties.
  
  \b Properties
  @li None
  
  \b Facilities
  @li \b eq_srcs Kinematic earthquake sources information.
  @li \b output Output manager associated with fault data.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from SingleRupture import SingleRupture
  eqsrcs = pyre.inventory.facilityArray("eq_srcs", itemFactory=eqsrcFactory,
                                        factory=SingleRupture)
  eqsrcs.meta['tip'] = "Kinematic earthquake sources information."
  
  from pylith.meshio.OutputFaultKin import OutputFaultKin
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputFaultKin)
  output.meta['tip'] = "Output manager associated with fault data."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoKi "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal_dir",
                     "final_slip",
                     "slip_time"],
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

    ModuleFaultCohesiveKin.quadrature(self, self.faultQuadrature)

    for eqsrc in self.eqsrcs.components():
      eqsrc.preinitialize()
    ModuleFaultCohesiveKin.eqsrcs(self, self.eqsrcs.inventory.facilityNames(),
                                  self.eqsrcs.components())

    for name in self.eqsrcs.inventory.facilityNames():
      self.availableFields['vertex']['info'] += ["final_slip_%s" % name]
      self.availableFields['vertex']['info'] += ["slip_time_%s" % name]

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
    ModuleFaultCohesiveKin.verifyConfiguration(self, self.mesh())

    for eqsrc in self.eqsrcs.components():
      eqsrc.verifyConfiguration()
    
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
    
    for eqsrc in self.eqsrcs.components():
      eqsrc.initialize()
    FaultCohesive.initialize(self, totalTime, numTimeSteps, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = ModuleFaultCohesiveKin.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveKin.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveKin.cellField(self, name)
    else:
      field = ModuleFaultCohesiveKin.cellField(self, name, fields)
    return field


  def finalize(self):
    """
    Cleanup.
    """
    for eqsrc in self.eqsrcs.components():
      eqsrc.finalize()
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
    self.eqsrcs = self.inventory.eqsrcs
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveKin.
    """
    ModuleFaultCohesiveKin.__init__(self)
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
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
