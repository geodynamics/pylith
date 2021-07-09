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
# @file pylith/apps/PetscApplication.py
#
# @brief Python PETSc application for creating an MPI application
# that uses PETSc.

from pythia.mpi import Application


class PetscApplication(Application):
    """Python PETSc application for creating an MPI application that uses PETSc.
    """

    import pythia.pyre.inventory

    # Dummy facility for passing options to PETSc
    from pylith.utils.PetscManager import PetscManager
    petsc = pythia.pyre.inventory.facility(
        "petsc", family="petsc_manager", factory=PetscManager)
    petsc.meta['tip'] = "Manager for PETSc options."

    includeCitations = pythia.pyre.inventory.bool(
        "include-citations", default=False)
    includeCitations.meta['tip'] = "At end of simulation, display information on how to cite PyLith and components used."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="petscapp"):
        """Constructor.
        """
        Application.__init__(self, name)
        return

    def onComputeNodes(self, *args, **kwds):
        """Run the application in parallel on the compute nodes.
        """
        self.petsc.initialize()

        if self.inventory.includeCitations:
            self.petsc.setOption("-citations", "")

            from pylith.utils.petsc import citationsRegister
            for entry in self.citations():
                citationsRegister(entry)

        try:
            self.main(*args, **kwds)
        except Exception as err:
            import traceback
            import sys

            self.cleanup()  # Attempt to clean up memory.
            print("Fatal error. Calling MPI_Abort() to abort PyLith application.")
            traceback.print_exc(file=sys.stdout)
            sys.stdout.flush()
            from pylith.mpi import mpi
            errorCode = -1
            mpi.mpi_abort(mpi.petsc_comm_world(), errorCode)

        self.cleanup()
        self.petsc.finalize()
        return

    def cleanup(self):
        """Deallocate data structures.
        """
        from pylith.utils.PetscComponent import PetscComponent
        for component in self.components():
            if isinstance(component, PetscComponent):
                component.cleanup()

            # Facility arrays are not PetscComponents but have components().
            elif hasattr(component, "components"):
                for subcomponent in component.components():
                    if isinstance(subcomponent, PetscComponent):
                        subcomponent.cleanup()

        self._cleanup()
        return

    def citations(self):
        """Register BibTeX entries for citing software.
        """
        return []

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        return

    def _cleanup(self):
        """Deallocate locally managed data structures.
        """
        return


# End of file
