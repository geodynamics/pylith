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

from .PropertyList import PropertyList
import pylith.utils.petsc as petsc


class PetscManager(PropertyList):
    """
    Manage PETSc options.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.petsc]
            ts_monitor = true
            ksp_monitor = true
            ksp_converged_reason = true
            snes_monitor = true
            snes_converged_reason = true
            snes_linesearch_monitor = true
        """
    }

    def __init__(self, name="petsc"):
        """Constructor.
        """
        PropertyList.__init__(self, name)

    def initialize(self):
        """Initialize PETSc.
        """
        import sys
        args = [sys.executable]
        options = self._getOptions()
        if len(options) > 0:
            for arg in options:
                args.append(arg)
        petsc.initialize(args)
        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()
        if 0 == comm.rank:
            self._info.log("Initialized PETSc.")

    def finalize(self):
        """Finalize PETSc.
        """
        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()
        if 0 == comm.rank:
            self._info.log("Finalizing PETSc.")
        petsc.finalize()

    def setOption(self, name, value):
        """Set option after PETSc initialization.
        """
        petsc.optionsSetValue(name, value)

    def _getOptions(self):
        """Cleanup options for PETSc.
        """
        args = []
        for iname, descriptor in self.items:
            args.append('-' + iname)
            if descriptor.value != 'true':
                args.append(descriptor.value)
        return args


# FACTORIES ////////////////////////////////////////////////////////////

def property_list():
    """
    Factory associated with PetscManager.
    """
    return PetscManager()


# End of file
