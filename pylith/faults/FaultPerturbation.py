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
#
# @file pylith/faults/TractionPerturbation.py
#
# @brief Python object for managing a traction perturbation.
#
# Factory: traction_perturbation

from pylith.utils.PetscComponent import PetscComponent
from .faults import TractionPerturbation as ModuleTractionPerturbation
from pylith.utils.NullComponent import NullComponent


def validateDir(value):
    """Validate direction.
    """
    msg = "Direction must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = list(map(float, value))
    except:
        raise ValueError(msg)
    return nums


class TractionPerturbation(PetscComponent, ModuleTractionPerturbation):
    """
    Time-dependent perturbation in fault traction. Implements `PetscComponent`.

    Traction perturbation can be initial fault tractions or a more complex spatial and temporal variation.

    :::{important}
    The components are specified in the local normal-tangential coordinate system for the fault.
    Ambiguities in specifying the shear (tangential) tractions in 3D problems are resolved using the `ref_dir_1` and `ref_dir_2` properties.
    The first tangential direction is $\\vec{z} \\times \\vec{r}_1$ unless these are colinear, then $\\vec{r}_2$ (`ref_dir_2`) is used.
    The second tangential direction is $\\vec{n} \\times \\vec{t}_1$.
    :::

    :::{seealso}
    See [`AuxSubfieldsTimeDependent` Component](AuxSubfieldsTimeDependent.md) for the functional form of the time dependence.
    :::
    """
    DOC_CONFIG = {
        "cfg": """
            # Perturbation in fault traction.
            [pylithapp.problem.interfaces.fault.traction_perturbations.perturbation]

            use_initial = False
            use_time_history = True
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Fault traction perturbation
            db_auxiliary_field.values = [time_history_amplitude_tangential, time_history_amplitude_normal, time_history_start_time]
            db_auxiliary_field.data = [2.0*MPa, -1.0*MPa, 0.0]

            time_history = spatialdata.spatialdb.TimeHistory
            time_history.description = Impulse time history
            time_history.filename = impulse.timedb
            """,
    }


    import pythia.pyre.inventory

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    def __init__(self, name="tractionperturbation"):
        """Constructor.
        """
        PetscComponent.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
        self.auxiliarySubfields = AuxSubfieldsTimeDependent("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log(
                "Performing minimal initialization of fault traction perturbation '%s'." % self.aliases[-1])

        PetscComponent.preinitialize(self, problem)

        ModuleTractionPerturbation.useInitial(self, self.useInitial)
        ModuleTractionPerturbation.useRate(self, self.useRate)
        ModuleTractionPerturbation.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleTractionPerturbation.setTimeHistoryDB(self, self.dbTimeHistory)
        return

    def _validate(self, context):
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            trait = self.inventory.getTrait("time_history")
            self._validationError(context, trait,
                f"Missing time history database for time-dependent fault traction perturbation '{self.aliases[-1]}'.")
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log(
                f"Time history for time-dependent fault traction perturbation '{self.aliases[-1]}' not enabled. Ignoring provided time history database.")

    def _validationError(self, context, trait, msg):
        from pythia.pyre.inventory.Item import Item
        error = ValueError(msg)
        descriptor = self.getTraitDescriptor(trait.name)
        context.error(error, items=[Item(trait, descriptor)])

    def _configure(self):
        """Setup members using inventory.
        """
        PetscComponent._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleTractionPerturbation.__init__(self)
        return


# Factories

def traction_perturbation():
    """Factory associated with TractionPerturbation.
    """
    return TractionPerturbation()


# End of file
