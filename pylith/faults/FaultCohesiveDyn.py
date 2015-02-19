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

## @file pylith/faults/FaultCohesiveDyn.py
##

## @brief Python object for a fault surface with dynamic
## (friction) fault implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveDyn as ModuleFaultCohesiveDyn

from pylith.utils.NullComponent import NullComponent

# FaultCohesiveDyn class
class FaultCohesiveDyn(FaultCohesive, Integrator, ModuleFaultCohesiveDyn):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesiveDyn facilities and properties.
  
  \b Properties
  @li \b zero_tolerance Tolerance for detecting zero values.
  @li \b open_free_surface If True, enforce traction free surface when
    the fault opens, otherwise use initial tractions even when the
    fault opens.
  
  \b Facilities
  @li \b tract_perturbation Prescribed perturbation in fault tractions.
  @li \b friction Fault constitutive model.
  @li \b output Output manager associated with fault data.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  zeroTolerance = pyre.inventory.float("zero_tolerance", default=1.0e-10,
                                       validator=pyre.inventory.greaterEqual(0.0))
  zeroTolerance.meta['tip'] = "Tolerance for detecting zero values."

  openFreeSurf = pyre.inventory.bool("open_free_surface", default=True)
  openFreeSurf.meta['tip'] = "If True, enforce traction free surface when " \
    "the fault opens, otherwise use initial tractions even when the " \
    "fault opens."

  tract = pyre.inventory.facility("traction_perturbation", family="traction_perturbation",
                               factory=NullComponent)
  tract.meta['tip'] = "Prescribed perturbation in fault tractions."

  from pylith.friction.StaticFriction import StaticFriction
  friction = pyre.inventory.facility("friction", family="friction_model",
                                     factory=StaticFriction)

  from pylith.meshio.OutputFaultDyn import OutputFaultDyn
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputFaultDyn)
  output.meta['tip'] = "Output manager associated with fault data."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivedyn"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoDy "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal_dir"],
            'data': ["slip",
                     "slip_rate",
                     "traction"]},
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

    ModuleFaultCohesiveDyn.quadrature(self, self.faultQuadrature)

    if mesh.dimension() == 2:
      self.availableFields['vertex']['info'] += ["strike_dir"]
    elif mesh.dimension() == 3:
      self.availableFields['vertex']['info'] += ["strike_dir",
                                                 "dip_dir"]

    if not isinstance(self.tract, NullComponent):
      self.tract.preinitialize(mesh)
      self.availableFields['vertex']['info'] += self.tract.availableFields['vertex']['info']

    self.availableFields['vertex']['info'] += \
        self.friction.availableFields['vertex']['info']
    self.availableFields['vertex']['data'] += \
        self.friction.availableFields['vertex']['data']
    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    FaultCohesive.verifyConfiguration(self)
    Integrator.verifyConfiguration(self)
    ModuleFaultCohesiveDyn.verifyConfiguration(self, self.mesh())

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
      field = ModuleFaultCohesiveDyn.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveDyn.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveDyn.cellField(self, name)
    else:
      field = ModuleFaultCohesiveDyn.cellField(self, name, fields)
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
    if not isinstance(self.inventory.tract, NullComponent):
      ModuleFaultCohesiveDyn.tractPerturbation(self, self.inventory.tract)
    ModuleFaultCohesiveDyn.frictionModel(self, self.inventory.friction)
    ModuleFaultCohesiveDyn.zeroTolerance(self, self.inventory.zeroTolerance)
    ModuleFaultCohesiveDyn.openFreeSurf(self, self.inventory.openFreeSurf)
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveDyn.
    """
    ModuleFaultCohesiveDyn.__init__(self)
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
  Factory associated with FaultCohesiveDyn.
  """
  return FaultCohesiveDyn()


# End of file 
