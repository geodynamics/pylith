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
# @file pylith/problems/Problem.py
#
# @brief Python abstract base class for crustal dynamics problems.
#
# Factory: problem.

from pylith.utils.PetscComponent import PetscComponent
from .problems import Problem as ModuleProblem
from .problems import Physics
from pylith.utils.NullComponent import NullComponent
from .ProblemDefaults import ProblemDefaults


# Factories for items in facility arrays

def materialFactory(name):
    """
    Factory for material items.
    """
    from pyre.inventory import facility
    from pylith.materials.Elasticity import Elasticity
    return facility(name, family="material", factory=Elasticity)


def bcFactory(name):
    """
    Factory for boundary condition items.
    """
    from pyre.inventory import facility
    from pylith.bc.DirichletTimeDependent import DirichletTimeDependent
    return facility(name, family="boundary_condition", factory=DirichletTimeDependent)


def faultFactory(name):
    """
    Factory for fault items.
    """
    from pyre.inventory import facility
    from pylith.faults.FaultCohesiveKin import FaultCohesiveKin
    return facility(name, family="fault", factory=FaultCohesiveKin)


def observerFactory(name):
    """
    Factory for output items.
    """
    from pyre.inventory import facility
    from pylith.meshio.OutputSolnDomain import OutputSolnDomain
    return facility(name, family="observer", factory=OutputSolnDomain)


class Problem(PetscComponent, ModuleProblem):
    """
    Python abstract base class for crustal dynamics problems.

    INVENTORY

    Properties
      - *dimension* Spatial dimension of problem space.
      - *formulation* Formulation for equations ('quasistatic' or 'dynamic').
      - *solver* Type of solver to use.

    Facilities
      - *solution* Solution field.
      - *normalizer* Nondimensionalizer for problem.
      - *materials* Array of materials (governing equations) in the problem.
      - *bc* Array of boundary conditions.
      - *interfaces* Array of interior surfaces with relative displacement constraints or constitutive models.
      - *solution_observers* Array of observers for solution.
      - *gravity_field* Gravity field for problem (SpatialDB).
      - *defaults* Default options for problem.
    """

    import pyre.inventory
    from pylith.utils.EmptyBin import EmptyBin

    defaults = pyre.inventory.facility("defaults", family="problem_defaults", factory=ProblemDefaults)
    defaults.meta['tip'] = "Default options for problem."

    formulation = pyre.inventory.str("formulation", default="quasistatic",
                                     validator=pyre.inventory.choice(["quasistatic", "dynamic", "dynamic_imex"]))
    formulation.meta['tip'] = "Formulation for equations."

    solverChoice = pyre.inventory.str("solver", default="linear",
                                      validator=pyre.inventory.choice(["linear", "nonlinear"]))
    solverChoice.meta['tip'] = "Type of solver to use ['linear', 'nonlinear']."

    from .Solution import Solution
    solution = pyre.inventory.facility("solution", family="solution", factory=Solution)
    solution.meta['tip'] = "Solution field for problem."

    from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic
    normalizer = pyre.inventory.facility("normalizer", family="nondimensional", factory=NondimElasticQuasistatic)
    normalizer.meta['tip'] = "Nondimensionalizer for problem."

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facilityArray("materials", itemFactory=materialFactory, factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    bc = pyre.inventory.facilityArray("bc", itemFactory=bcFactory, factory=EmptyBin)
    bc.meta['tip'] = "Boundary conditions."

    interfaces = pyre.inventory.facilityArray("interfaces", itemFactory=faultFactory, factory=EmptyBin)
    interfaces.meta['tip'] = "Interior surfaces with constraints or constitutive models."

    from pylith.problems.SingleObserver import SingleSolnObserver
    observers = pyre.inventory.facilityArray(
        "solution_observers", itemFactory=observerFactory, factory=SingleSolnObserver)
    observers.meta['tip'] = "Observers (e.g., output) for solution."

    gravityField = pyre.inventory.facility("gravity_field", family="spatial_database", factory=NullComponent)
    gravityField.meta['tip'] = "Database used for gravity field."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="problem"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="problem")
        self.mesh = None
        return

    def preinitialize(self, mesh):
        """
        Do minimal initialization.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization before verifying configuration.")

        self._createModuleObj()
        ModuleProblem.setIdentifier(self, self.aliases[-1])
        self.defaults.preinitialize()

        if self.formulation == "quasistatic":
            formulationType = Physics.QUASISTATIC
        elif self.formulation == "dynamic":
            formulationType = Physics.DYNAMIC
        elif self.formulation == "dynamic_imex":
            formulationType = Physics.DYNAMIC_IMEX
        else:
            raise ValueError("Unknown formulation '{}'.".format(self.formulation))
        ModuleProblem.setFormulation(self, formulationType)

        if self.solverChoice == "linear":
            ModuleProblem.setSolverType(self, ModuleProblem.LINEAR)
        elif self.solverChoice == "nonlinear":
            ModuleProblem.setSolverType(self, ModuleProblem.NONLINEAR)
        else:
            raise ValueError("Unknown solver choice '%s'." % self.solverChoice)
        ModuleProblem.setNormalizer(self, self.normalizer)
        if not isinstance(self.gravityField, NullComponent):
            ModuleProblem.setGravityField(self, self.gravityField)

        # Do minimal setup of solution.
        self.solution.preinitialize(self, mesh)
        ModuleProblem.setSolution(self, self.solution.field)

        # Preinitialize materials
        for material in self.materials.components():
            material.preinitialize(self)
        ModuleProblem.setMaterials(self, self.materials.components())

        # Preinitialize boundary conditions.
        for bc in self.bc.components():
            bc.preinitialize(self)
        ModuleProblem.setBoundaryConditions(self, self.bc.components())

        # Preinitialize interfaces
        for interface in self.interfaces.components():
            interface.preinitialize(self)
        ModuleProblem.setInterfaces(self, self.interfaces.components())

        # Preinitialize observers.
        for observer in self.observers.components():
            observer.preinitialize(self)
            ModuleProblem.registerObserver(self, observer)

        ModuleProblem.preinitialize(self, mesh)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Verifying compatibility of problem configuration.")

        ModuleProblem.verifyConfiguration(self)

        self._printInfo()
        return

    def initialize(self):
        """
        Initialize integrators and constraints.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Initializing {} problem.".format(self.formulation))

        ModuleProblem.initialize(self)
        return

    def run(self, app):
        """
        Solve the problem.
        """
        raise NotImplementedError("run() not implemented.")
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Finalizing problem.")
        return

    def checkpoint(self):
        """
        Save problem state for restart.
        """
        raise NotImplementedError("checkpoint() not implemented.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _printInfo(self):
        """Write overview of problem to info journal.
        """
        msg = (
            "Scales for nondimensionalization:",
            "    Length scale: {}".format(self.normalizer.getLengthScale()),
            "    Time scale: {}".format(self.normalizer.getTimeScale()),
            "    Pressure scale: {}".format(self.normalizer.getPressureScale()),
            "    Density scale: {}".format(self.normalizer.getDensityScale()),
            "    Temperature scale: {}".format(self.normalizer.getTemperatureScale()),
        )
        self._info.log("\n".join(msg))
        return

    def _setupLogging(self):
        """
        Setup event logging.
        """
        if not "_loggingPrefix" in dir(self):
            self._loggingPrefix = ""

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("Problem")
        logger.initialize()

        self._eventLogger = logger
        return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """
    Factory associated with Problem.
    """
    return Problem()


# End of file
