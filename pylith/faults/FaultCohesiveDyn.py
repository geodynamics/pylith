# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveDyn as ModuleFaultCohesiveDyn


def perturbationFactory(name):
    """
    Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .TractionPerturbation import TractionPerturbation
    return facility(name, family="traction_perturbation", factory=TractionPerturbation)


class FaultCohesiveDyn(FaultCohesive, ModuleFaultCohesiveDyn):
    """
    Fault surface with dynamic (spontaneous) slip implemented with cohesive cells.

    The fault may have an arbitrary number of traction perturbations.
    They are superimposed at each time step to create prescribed tractions on the fault.

    Implements `FaultCohesive`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Specify spontaneous rupture on a fault with static friction and a traction perturbation.
            [pylithapp.problem.interfaces.fault]
            label = fault
            edge = fault_edge

            observers.observer.data_fields = [slip]

            traction_perturbations = [initial]
            fault_rheology = pylith.faults.FrictionStatic

            [pylithapp.problem.interfaces.fault.traction_perturbations.initial]
            origin_time = 1*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Fault rupture auxiliary field spatial database
            db_auxiliary_field.values = [traction_left_lateral, traction_opening]
            db_auxiliary_field.data = [-2.0*MPa, -10.0*MPa]

            [pylithapp.problem.interfaces.fault.fault_rheology]
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Static friction spatial database
            db_auxiliary_field.values = [cohesion, static_coefficient]
            db_auxiliary_field.data = [0.0*Pa, 0.5]
            """
    }

    import pythia.pyre.inventory

    from pylith.utils.EmptyBin import EmptyBin
    tractionPerturbations = pythia.pyre.inventory.facilityArray("traction_perturbations", itemFactory=perturbationFactory, factory=EmptyBin)
    tractionPerturbations.meta['tip'] = "Fault traction perfurbations."

    from pylith.faults.FrictionStatic import FrictionStatic
    faultRheology = pythia.pyre.inventory.facility("fault_rheology", family="fault_rheology", factory=FrictionStatic)
    faultRheology.meta['tip'] = "Fault constitutive model."

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    def __init__(self, name="faultcohesivedyn"):
        """Initialize configuration.
        """
        FaultCohesive.__init__(self, name)

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log("Pre-initializing fault '%s'." % self.labelName)

        FaultCohesive.preinitialize(self, problem)

        for perturbation in self.tractionPerturbations.components():
            perturbation.preinitialize(problem)
        ModuleFaultCohesiveDyn.setTractionPerturbations(
            self, self.tractionPerturbations.inventory.facilityNames(), self.tractionPerturbations.components())

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveDyn.verifyConfiguration(self, self.mesh())

        for perturbation in self.tractionPerturbations.components():
            perturbation.verifyConfiguration()

    def finalize(self):
        """Cleanup.
        """
        for perturbation in self.tractionPerturbations.components():
            perturbation.finalize()
        FaultCohesive.finalize(self)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        FaultCohesive._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesiveDyn.
        """
        ModuleFaultCohesiveDyn.__init__(self)
        return


# Factories

def fault():
    """Factory associated with FaultCohesiveDyn.
    """
    return FaultCohesiveDyn()


# End of file
