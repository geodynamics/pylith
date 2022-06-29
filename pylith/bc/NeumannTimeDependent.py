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
# @file pylith/bc/NeumannTimeDependent.py
#
# @brief Python object for managing a time-dependent Neumann (natural) boundary condition.
#
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from .bc import NeumannTimeDependent as ModuleNeumannTimeDependent
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


class NeumannTimeDependent(BoundaryCondition, ModuleNeumannTimeDependent):
    """
    Neumann time-dependent boundary condition. Implements `BoundaryCondition`.

    This boundary condition applies a Neumann boundary condition for a single solution subfield on a boundary.
    To apply Neumann boundary conditions for multiple solution subfields on a boundary, use multiple Neumann boundary conditions.

    :::{important}
    The components are specified in the local normal-tangential coordinate system for the boundary. Ambiguities in specifying the shear (tangential) tractions in 3D problems are resolved using the `ref_dir_1` and `ref_dir_2` properties.
    The first tangential direction is $\\vec{z} \\times \\vec{r}_1$ unless these are colinear, then $\\vec{r}_2$ (`ref_dir_2`) is used.
    The second tangential direction is $\\vec{n} \\times \\vec{t}_1$.
    :::

    :::{seealso}
    See [`AuxSubfieldsTimeDependent` Component](AuxSubfieldsTimeDependent.md) for the functional form of the time depenence.
    :::
    """
    DOC_CONFIG = {
        "cfg": """
            # Neumann (traction) boundary condition in 2D on -y boundary.
            [pylithapp.problem.bc.bc_yneg]
            label = boundary_yneg
            field = displacement
            scale_name = pressure

            use_initial = False
            use_time_history = True
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Displacement Neumann BC +y boundary
            db_auxiliary_field.values = [time_history_amplitude_tangential, time_history_amplitude_normal, time_history_start_time]
            db_auxiliary_field.data = [2.0*MPa, -1.0*MPa, 0.0]

            time_history = spatialdata.spatialdb.TimeHistory
            time_history.description = Impulse time history
            time_history.filename = impulse.timedb
            """,
    }


    import pythia.pyre.inventory

    scaleName = pythia.pyre.inventory.str("scale_name", default="pressure",
                                   validator=pythia.pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
    scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    refDir1 = pythia.pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3D."

    refDir2 = pythia.pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3D."

    def __init__(self, name="neumanntimedependent"):
        """Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
        self.auxiliarySubfields = AuxSubfieldsTimeDependent("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log(
                "Performing minimal initialization of time-dependent Neumann boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition.preinitialize(self, problem)

        ModuleNeumannTimeDependent.setScaleName(self, self.scaleName)
        ModuleNeumannTimeDependent.setRefDir1(self, self.refDir1)
        ModuleNeumannTimeDependent.setRefDir2(self, self.refDir2)
        ModuleNeumannTimeDependent.useInitial(self, self.useInitial)
        ModuleNeumannTimeDependent.useRate(self, self.useRate)
        ModuleNeumannTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleNeumannTimeDependent.setTimeHistoryDB(self, self.dbTimeHistory)
        return

    def _validate(self, context):
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            trait = self.inventory.getTrait("time_history")
            self._validationError(context, trait,
                f"Missing time history database for time-dependent Neumann boundary condition '{self.aliases[-1]}'.")
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log(
                f"Time history for time-dependent Neumann boundary condition '{self.aliases[-1]}' not enabled. Ignoring provided time history database.")

    def _validationError(self, context, trait, msg):
        from pythia.pyre.inventory.Item import Item
        error = ValueError(msg)
        descriptor = self.getTraitDescriptor(trait.name)
        context.error(error, items=[Item(trait, descriptor)])

    def _configure(self):
        """Setup members using inventory.
        """
        BoundaryCondition._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleNeumannTimeDependent.__init__(self)
        return


# Factories

def boundary_condition():
    """Factory associated with NeumannTimeDependent.
    """
    return NeumannTimeDependent()


# End of file
