# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/utils/PetscManager.py
#
# @brief Python PetscManager object for managing PETSc options.
#
# The PetscManager also takes care of initializing and finalizing
# PETSc.
#
# Factory: petsc_manager

from pythia.pyre.components.Component import Component
import pylith.utils.petsc as petsc


class PetscManager(Component):
    """
    Python PetscManager object for managing PETSc options.

    Factory: petsc_manager
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="petsc"):
        """
        Constructor.
        """
        Component.__init__(self, name, facility="petsc_manager")
        self.options = []
        return

    def initialize(self):
        """
        Initialize PETSc.
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
        return

    def finalize(self):
        """
        Finalize PETSc.
        """
        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()
        if 0 == comm.rank:
            self._info.log("Finalizing PETSc.")
        petsc.finalize()
        return

    def setOption(self, name, value):
        """
        Set option after PETSc initialization.
        """
        petsc.optionsSetValue(name, value)
        return

    def updateConfiguration(self, registry):
        """
        Update Pyre configuration.
        """
        self.options = [
            (name, descriptor) for name, descriptor in registry.properties.items()
        ]
        return []

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _getOptions(self):
        """
        Cleanup options for PETSc.
        """
        args = []
        for iname, descriptor in self.options:
            args.append('-' + iname)
            if descriptor.value != 'true':
                args.append(descriptor.value)
        return args


# FACTORIES ////////////////////////////////////////////////////////////

def petsc_manager():
    """
    Factory associated with PetscManager.
    """
    return PetscManager()


# End of file
