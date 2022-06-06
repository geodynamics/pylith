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

from pylith.utils.PetscComponent import PetscComponent
from .problems import InitialCondition as ModuleInitialCondition


class InitialCondition(PetscComponent, ModuleInitialCondition):
    """
    Abstract base class for specifying initial conditions for the solution.
    """
    import pythia.pyre.inventory

    subfields = pythia.pyre.inventory.list("subfields", default=["displacement"])
    subfields.meta["tip"] = "Names of solution subfields for initial condition."

    def __init__(self, name="initialconditions"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="initial_conditions")

    def preinitialize(self, problem):
        """Setup initial conditions.
        """
        self._createModuleObj()
        ModuleInitialCondition.setIdentifier(self, self.aliases[-1])
        ModuleInitialCondition.setSubfields(self, self.subfields)

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")


# End of file
