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
#
# @file pylith/utils/PetscManager.py
#
# @brief Python PetscManager object for managing PETSc options.
#
# The PetscManager also takes care of initializing and finalizing
# PETSc.
#
# Factory: petsc_manager

from .PropertyList import PropertyList
import pylith.utils.petsc as petsc


class PetscManager(PropertyList):
    """Python PetscManager object for managing PETSc options.

    Factory: property_list
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="petsc"):
        """Constructor.
        """
        PropertyList.__init__(self, name)
        return

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
        return

    def finalize(self):
        """Finalize PETSc.
        """
        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()
        if 0 == comm.rank:
            self._info.log("Finalizing PETSc.")
        petsc.finalize()
        return

    def setOption(self, name, value):
        """Set option after PETSc initialization.
        """
        petsc.optionsSetValue(name, value)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

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
