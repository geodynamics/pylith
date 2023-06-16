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
            monitors = True
            parallel = False
            initial_guess = True
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

    testing = pythia.pyre.inventory.bool("testing", default=False)
    testing.meta["tip"] = "Use default PETSc testing options."

    def __init__(self, name="petscdefaults"):
        """Constructor.
        """
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
        if self.testing:
            value |= ModuleDefaults.TESTING
        return value


# FACTORIES ////////////////////////////////////////////////////////////

def petsc_defaults():
    """Factory associated with PetscDefaults.
    """
    return PetscDefaults()


# End of file
