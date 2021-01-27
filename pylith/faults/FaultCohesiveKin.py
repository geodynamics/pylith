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
from .faults import FaultCohesiveKin as ModuleFaultCohesiveKin

# ITEM FACTORIES ///////////////////////////////////////////////////////


def eqsrcFactory(name):
    """
    Factory for earthquake source items.
    """
    from pyre.inventory import facility
    from .KinSrcStep import KinSrcStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcStep)


class FaultCohesiveKin(FaultCohesive, ModuleFaultCohesiveKin):
    """
    Python object for a fault surface with kinematic (prescribed) slip
    implemented with cohesive elements.

    INVENTORY

    Properties
      - None

    Facilities
      - *eq_srcs* Kinematic earthquake sources information.
      - *observers* Observers of the fault (e.g., output).

    FACTORY: fault
    """

    import pyre.inventory

    from .SingleRupture import SingleRupture
    eqRuptures = pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupture)
    eqRuptures.meta['tip'] = "Kinematic earthquake sources information."

    from pylith.utils.EmptyBin import EmptyBin
    from pylith.topology.Subfield import subfieldFactory
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxiliarySubfields.meta['tip'] = "Discretization of fault parameters."

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    #from pylith.meshio.OutputFaultKin import OutputFaultKin
    #outputManager = pyre.inventory.facility("output", family="output_manager", factory=OutputFaultKin)
    #output.meta['tip'] = "Output manager associated with fault information."

    def __init__(self, name="faultcohesivekin"):
        """
        Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Pre-initializing fault '%s'." % self.label)

        FaultCohesive.preinitialize(self, problem)

        for eqsrc in self.eqRuptures.components():
            eqsrc.preinitialize()
        ModuleFaultCohesiveKin.setEqRuptures(
            self, self.eqRuptures.inventory.facilityNames(), self.eqRuptures.components())

        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveKin.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqRuptures.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """
        Cleanup.
        """
        for eqsrc in self.eqRuptures.components():
            eqsrc.finalize()
        FaultCohesive.finalize(self)
        # self.output.close()
        # self.output.finalize()
        return

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


# Factories

def fault():
    """
    Factory associated with FaultCohesiveKin.
    """
    return FaultCohesiveKin()


# End of file
