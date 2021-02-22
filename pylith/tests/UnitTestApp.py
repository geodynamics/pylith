# ======================================================================
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
# ======================================================================
#
# @file pylith/tests/UnitTestApp.py
#
# @brief Python application for Python unit tests.

from pythia.pyre.applications.Script import Script

import unittest


class UnitTestApp(Script):
    """
    Test application.
    """
    cov = None
    try:
        import coverage
        cov = coverage.Coverage(source=["spatialdata"])
    except ImportError:
        pass

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="unittestapp", petsc_options=[("malloc_dump", "true")]):
        """
        Constructor.
        """
        Script.__init__(self, name)
        self.petscOptions = petsc_options
        return

    def main(self):
        """
        Run the application.
        """
        if self.cov:
            self.cov.start()

        from pylith.utils.PetscManager import PetscManager
        petsc = PetscManager()
        petsc.options = self.petscOptions
        petsc.initialize()

        success = unittest.TextTestRunner(verbosity=2).run(self._suite()).wasSuccessful()

        petsc.finalize()

        if self.cov:
            self.cov.stop()
            self.cov.save()

        if not success:
            import sys
            sys.exit(1)
        return


def configureSubcomponents(facility):
    """Configure subcomponents."""
    for component in facility.components():
        configureSubcomponents(component)
        component._configure()
    return


# End of file
