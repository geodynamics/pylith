#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

"""Script to run pylith Python test suite.

Run `coverage report` to generate a report (included).
Run `coverage html -d DIR` to generate an HTML report in directory `DIR`.
"""


import unittest
import sys
import importlib


class TestApp:
    """Application to run tests.
    """
    DIRS = [
        "apps",
        "bc",
        "faults",
        "materials",
        "meshio",
        "mpi",
        "problems",
        "topology",
        "utils",
    ]

    def __init__(self):
        self.cov = None
        try:
            import coverage
            self.cov = coverage.Coverage(source=["pylith"])
        except ImportError:
            pass

    def run(self):
        """Run the application.
        """
        if self.cov:
            self.cov.start()

        from pylith.utils.PetscManager import PetscManager
        petsc = PetscManager()
        petsc.initialize()

        success = unittest.TextTestRunner(verbosity=2).run(self._suite()).wasSuccessful()

        petsc.finalize()

        if not success:
            sys.exit(1)

        if self.cov:
            self.cov.stop()
            self.cov.save()
            self.cov.report()
            self.cov.xml_report(outfile="coverage.xml")

    def _suite(self):
        """Setup the test suite.
        """
        for d in self.DIRS:
            sys.path.append(f"./{d}")

        test_modules = []
        for d in self.DIRS:
            mod = importlib.import_module(d)
            test_modules += mod.test_modules()

        loader = unittest.defaultTestLoader
        suite = unittest.TestSuite()
        for mod in test_modules:
            suite.addTests(loader.loadTestsFromModule(mod))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestApp().run()


# End of file
