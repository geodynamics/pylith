# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveKin as ModuleFaultCohesiveKin


def eqsrcFactory(name):
    """
    Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .KinSrcStep import KinSrcStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcStep)


class FaultCohesiveKin(FaultCohesive, ModuleFaultCohesiveKin):
    """
    Fault surface with kinematic (prescribed) slip implemented with cohesive cells.

    The fault may have an arbitrary number of kinematic sources for coseismic slip and creep.
    They are superimposed at each time step to create the prescribed slip on the fault.

    Implements `FaultCohesive`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Specify prescribed slip on a fault via two earthquakes in a 2D domain.
            [pylithapp.problem.interfaces.fault]
            label = fault
            id = 10

            observers.observer.data_fields = [slip]

            # Two earthquakes with different slip time functions.
            eq_ruptures = [quake10, quake50]
            quake10 = pylith.faults.KinSrcBrune
            quake50 = pylith.faults.KinSrcLiuCosine

            # Rupture parameters for the first earthquake.
            [pylithapp.problem.interfaces.fault.eq_ruptures.quake10]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.label = Fault rupture auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

            # Rupture parameters for the second earthquake.
            [pylithapp.problem.interfaces.fault.eq_ruptures.quake50]
            origin_time = 50*year
            
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.label = Fault rupture auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -1.0*m, 0.0*m]
            """
    }

    import pythia.pyre.inventory

    from .SingleRupture import SingleRupture
    eqRuptures = pythia.pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupture)
    eqRuptures.meta['tip'] = "Kinematic earthquake sources information."

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    #from pylith.meshio.OutputFaultKin import OutputFaultKin
    #outputManager = pythia.pyre.inventory.facility("output", family="output_manager", factory=OutputFaultKin)
    #output.meta['tip'] = "Output manager associated with fault information."

    def __init__(self, name="faultcohesivekin"):
        """Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
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
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveKin.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqRuptures.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """Cleanup.
        """
        for eqsrc in self.eqRuptures.components():
            eqsrc.finalize()
        FaultCohesive.finalize(self)
        # self.output.close()
        # self.output.finalize()
        return

    def _configure(self):
        """Setup members using inventory.
        """
        FaultCohesive._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesiveKin.
        """
        ModuleFaultCohesiveKin.__init__(self)
        return


# Factories

def fault():
    """Factory associated with FaultCohesiveKin.
    """
    return FaultCohesiveKin()


# End of file
