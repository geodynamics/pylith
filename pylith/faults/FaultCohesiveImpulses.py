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

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveImpulses as ModuleFaultCohesiveImpulses

class FaultCohesiveImpulses(FaultCohesive, ModuleFaultCohesiveImpulses):
    """Python object for a fault surface with slip impulses for Green's
    functions implemented with cohesive elements.

    Factory: fault
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pythia.pyre.inventory
    from pythia.pyre.units.length import m

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    threshold = pythia.pyre.inventory.dimensional("threshold", default=1.0e-6*m)
    threshold.meta['tip'] = "Threshold for non-zero amplitude."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="faultcohesiveimpulses"):
        """
        Initialize configuration.
        """
        FaultCohesive.__init__(self, name)


    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log("Pre-initializing fault '%s'." % self.label())
        FaultCohesive.preinitialize(self, problem)
        #ModuleFaultCohesiveImpulses.setThreshold(self.threshold.value)
  

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveImpulses.verifyConfiguration(self, self.mesh())

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesiveImpulses.
        """
        ModuleFaultCohesiveImpulses.__init__(self)
  

# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveImpulses.
  """
  return FaultCohesiveImpulses()


# End of file 
