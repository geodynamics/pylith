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

from .BoundaryCondition import BoundaryCondition
from .bc import DirichletTimeDependent as ModuleDirichletTimeDependent
from pylith.utils.NullComponent import NullComponent


class DirichletTimeDependent(BoundaryCondition, ModuleDirichletTimeDependent):
    """
    Dirichlet (prescribed values) time-dependent boundary condition.

    This boundary condition sets values of a single solution subfield on a boundary.
    To set multiple solution subfields on a boundary, use multiple Dirichlet boundary conditions.

    :::{seealso}
    See [`AuxSubfieldsTimeDependent` Component](AuxSubfieldsTimeDependent.md) for the functional form of the time depenence.
    :::

    Implements `BoundaryCondition`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Dirichlet (prescribed displacements) boundary condition constraining the x and y degrees of freedom on the +y boundary.
            [pylithapp.problem.bc.bc_ypos]
            constrained_dof = [0, 1]
            label = boundary_ypos
            field = displacement

            use_initial = False
            use_time_history = True
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Displacement Dirichlet BC +y boundary
            db_auxiliary_field.values = [time_history_amplitude_x, time_history_amplitude_y, time_history_start_time]
            db_auxiliary_field.data = [1.0*m, 0.0*m, 0.0]

            time_history = spatialdata.spatialdb.TimeHistory
            time_history.description = Impulse time history
            time_history.filename = impulse.timedb
            """,
    }


    import pythia.pyre.inventory

    constrainedDOF = pythia.pyre.inventory.array("constrained_dof", converter=int, default=[])
    constrainedDOF.meta['tip'] = "Array of constrained degrees of freedom (0=1st DOF, 1=2nd DOF, etc)."

    useInitial = pythia.pyre.inventory.bool("use_initial", default=True)
    useInitial.meta['tip'] = "Use initial term in time-dependent expression."

    useRate = pythia.pyre.inventory.bool("use_rate", default=False)
    useRate.meta['tip'] = "Use rate term in time-dependent expression."

    useTimeHistory = pythia.pyre.inventory.bool("use_time_history", default=False)
    useTimeHistory.meta['tip'] = "Use time history term in time-dependent expression."

    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=NullComponent, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude."

    def __init__(self, name="dirichlettimedependent"):
        """Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsTimeDependent import AuxSubfieldsTimeDependent
        self.auxiliarySubfields = AuxSubfieldsTimeDependent(
            "auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        import numpy

        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log(
                "Performing minimal initialization of time-dependent Dirichlet boundary condition '%s'." % self.aliases[-1])

        BoundaryCondition.preinitialize(self, problem)

        ModuleDirichletTimeDependent.setConstrainedDOF(
            self, numpy.array(self.constrainedDOF, dtype=numpy.int32))
        ModuleDirichletTimeDependent.useInitial(self, self.useInitial)
        ModuleDirichletTimeDependent.useRate(self, self.useRate)
        ModuleDirichletTimeDependent.useTimeHistory(self, self.useTimeHistory)
        if not isinstance(self.dbTimeHistory, NullComponent):
            ModuleDirichletTimeDependent.setTimeHistoryDB(
                self, self.dbTimeHistory)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        BoundaryCondition.verifyConfiguration(self, self.mesh())
        spaceDim = self.mesh().coordsys().getSpaceDim()
        for d in self.bcDOF:
            if d < 0 or d >= spaceDim:
                raise ValueError("Attempting to constrain DOF (%d) that doesn't exist for time-dependent Dirichlet boundary condition '%s'. Space dimension is %d." %
                                 (d, self.aliases[-1], spaceDim))
        return

    def _validate(self, context):
        if 0 == len(self.constrainedDOF):
            trait = self.inventory.getTrait("constrained_dof")
            self._validationError(context, trait, f"No constrained degrees of freedom found for time-dependent Dirichlet boundary condition '{self.aliases[-1]}'. "
                "'constrained_dof' must be a zero-based integer array (0=x, 1=y, 2=z).")
        if self.inventory.useTimeHistory and isinstance(self.inventory.dbTimeHistory, NullComponent):
            trait = self.inventory.getTrait("time_history")
            self._validationError(context, trait,
                f"Missing time history database for time-dependent Dirichlet boundary condition '{self.aliases[-1]}'.")
        if not self.inventory.useTimeHistory and not isinstance(self.inventory.dbTimeHistory, NullComponent):
            self._warning.log(
                f"Time history for time-dependent Dirichlet boundary condition '{self.aliases[-1]}' not enabled. Ignoring provided time history database.")

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
        ModuleDirichletTimeDependent.__init__(self)
        return


# Factories

def boundary_condition():
    """Factory associated with DirichletTimeDependent.
    """
    return DirichletTimeDependent()


# End of file
