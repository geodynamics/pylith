# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component


class PetscDefaults(Component):
    """
    Flags controlling use of default PETSc settings.
    No user-specified settings will be overwritten.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.petsc_defaults]
            solver = True
            parallel = False
            monitors = True
            initial_guess = True
            adaptive_time_stepping = False
            collective_io = True
            testing = False
        """
    }

    import pythia.pyre.inventory

    solver = pythia.pyre.inventory.bool("solver", default=True)
    solver.meta["tip"] = "Use default solver settings based on governing equations."

    parallel = pythia.pyre.inventory.bool("parallel", default=False)
    parallel.meta["tip"] = "Use solver settings normally used when running in parallel."

    monitors = pythia.pyre.inventory.bool("monitors", default=True)
    monitors.meta["tip"] = "Use default solver monitors."

    initialGuess = pythia.pyre.inventory.bool("initial_guess", default=True)
    initialGuess.meta["tip"] = "Use initial guess options."

    collectiveIO = pythia.pyre.inventory.bool("collective_io", default=True)
    collectiveIO.meta["tip"] = "Use default PETSc collective I/O options."

    testing = pythia.pyre.inventory.bool("testing", default=False)
    testing.meta["tip"] = "Use default PETSc testing options."

    adaptiveTS = pythia.pyre.inventory.bool("adaptive_time_stepping", default=False)
    adaptiveTS.meta["tip"] = "Use adaptive time stepping."

    impulseTS = pythia.pyre.inventory.bool("impulse_ts", default=False)
    impulseTS.meta["tip"] = "Use custom PETSc TS for impulses."

    def __init__(self, name="petscdefaults"):
        """Constructor."""
        Component.__init__(self, name)

    def flags(self):
        from .utils import PetscDefaults as ModuleDefaults

        value = ModuleDefaults.NONE
        if self.solver:
            value |= ModuleDefaults.SOLVER
        if self.parallel:
            value |= ModuleDefaults.PARALLEL
        if self.monitors:
            value |= ModuleDefaults.MONITORS
        if self.initialGuess:
            value |= ModuleDefaults.INITIAL_GUESS
        if self.collectiveIO:
            value |= ModuleDefaults.COLLECTIVE_IO
        if self.testing:
            value |= ModuleDefaults.TESTING
        if self.adaptiveTS:
            value |= ModuleDefaults.TS_ADAPTIVE
        if self.impulseTS:
            value |= ModuleDefaults.TS_IMPULSE
        return value


# FACTORIES ////////////////////////////////////////////////////////////


def petsc_defaults():
    """Factory associated with PetscDefaults."""
    return PetscDefaults()


# End of file
