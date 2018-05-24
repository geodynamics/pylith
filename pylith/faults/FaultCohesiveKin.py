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

# @file pylith/faults/FaultCohesiveKin.py
#

# @brief Python object for a fault surface with kinematic
# (prescribed) slip implemented with cohesive elements.
#
# Factory: fault

from .FaultCohesive import FaultCohesive
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise
from .faults import FaultCohesiveKin as ModuleFaultCohesiveKin

# ITEM FACTORIES ///////////////////////////////////////////////////////


def eqsrcFactory(name):
    """
    Factory for earthquake source items.
    """
    from pyre.inventory import facility
    from .KinSrc import KinSrc
    return facility(name, family="eq_kinematic_src", factory=KinSrc)


# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive, IntegratorPointwise, ModuleFaultCohesiveKin):
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

    Factory: fault
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    from SingleRupture import SingleRupture
    eqsrcs = pyre.inventory.facilityArray("eq_srcs", itemFactory=eqsrcFactory, factory=SingleRupture)
    eqsrcs.meta['tip'] = "Kinematic earthquake sources information."

    #from pylith.meshio.OutputFaultKin import OutputFaultKin
    #outputManager = pyre.inventory.facility("output", family="output_manager", factory=OutputFaultKin)
    #output.meta['tip'] = "Output manager associated with fault information."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="faultcohesivekin"):
        """
        Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        IntegratorPointwise.__init__(self)
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
        Integratorpointwise.preinitialize(self, mesh)

        for eqsrc in self.eqsrcs.components():
            eqsrc.preinitialize()
        ModuleFaultCohesiveKin.eqsrcs(self, self.eqsrcs.inventory.facilityNames(), self.eqsrcs.components())

        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        Integrator.verifyConfiguration(self)
        ModuleFaultCohesiveKin.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqsrcs.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """
        Cleanup.
        """
        for eqsrc in self.eqsrcs.components():
            eqsrc.finalize()
        FaultCohesive.finalize(self)
        Integrator.finalize(self)
        # self.output.close()
        # self.output.finalize()
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
        Create handle to C++ FaultCohesiveKin.
        """
        ModuleFaultCohesiveKin.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def fault():
    """
    Factory associated with FaultCohesiveKin.
    """
    return FaultCohesiveKin()


# End of file
