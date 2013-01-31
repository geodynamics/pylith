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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/faults/FaultCohesiveDynKin.py
##

## @brief Python object for a fault surface with dynamic
## (friction) fault implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive
from pylith.feassemble.Integrator import Integrator
from faults import FaultCohesiveDynKin as ModuleFaultCohesiveDynKin

from pylith.utils.NullComponent import NullComponent

# ITEM FACTORIES //////////////////////////////////////////////////////////

def eqsrcFactory(name):
  """
  Factory for earthquake source items.
  """
  from pyre.inventory import facility
  from EqKinSrc import EqKinSrc
  return facility(name, family="eq_kinematic_src", factory=EqKinSrc)

# FaultCohesiveDynKin class
class FaultCohesiveDynKin(FaultCohesive, Integrator, ModuleFaultCohesiveDynKin):
  """
  Python object for a fault surface with kinematic(prescribed) slip
  implemented with cohesive elements and dynamic friction.

  Inventory

  @class Inventory
  Python object for managing FaultCohesiveDynKin facilities and properties.
  
  \b Properties
  @li \b zero_tolerance Tolerance for detecting zero values.
  @li \b open_free_surface If True, enforce traction free surface when
    the fault opens, otherwise use initial tractions even when the
    fault opens.
  
  \b Facilities
  @li \b tract_perturbation Prescribed perturbation in fault tractions.
  @li \b dkSelector Dynamic Kinematic Selector.
  @li \b friction Fault constitutive model.
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

  from pylith.meshio.OutputFaultDynKin import OutputFaultDynKin
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputFaultDynKin)
  output.meta['tip'] = "Output manager associated with fault data."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivedynkin"):
    """
    Initialize configuration.
    """
    FaultCohesive.__init__(self, name)
    Integrator.__init__(self)
    self._loggingPrefix = "CoDyKi "

    self.availableFields = \
        {'vertex': \
           {'info': ["normal_dir",
	 	     "final_slip",
		     "slip_time"],
            'data': ["slip",
                     "slip_rate",
                     "traction",
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

    ModuleFaultCohesiveDynKin.quadrature(self, self.faultQuadrature)

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
    ModuleFaultCohesiveDyn.verifyConfiguration(self, self.mesh)

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
      field = ModuleFaultCohesiveDynKin.vertexField(self, name)
    else:
      field = ModuleFaultCohesiveDynKin.vertexField(self, name, fields)
    return field


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    if None == fields:
      field = ModuleFaultCohesiveDynKin.cellField(self, name)
    else:
      field = ModuleFaultCohesiveDynKin.cellField(self, name, fields)
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
      ModuleFaultCohesiveDynKin.tractPerturbation(self, self.inventory.tract)
    ModuelFaultCohesiveDynKin.dkSelector(self, self.inventory.dksel)
    ModuleFaultCohesiveDynKin.frictionModel(self, self.inventory.friction)
    ModuleFaultCohesiveDynKin.zeroTolerance(self, self.inventory.zeroTolerance)
    ModuleFaultCohesiveDynKin.openFreeSurf(self, self.inventory.openFreeSurf)
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to C++ FaultCohesiveDyn.
    """
    ModuleFaultCohesiveDynKin.__init__(self)
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
  return FaultCohesiveDynKin()


# End of file 
