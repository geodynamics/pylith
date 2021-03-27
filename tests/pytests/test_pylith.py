#!/usr/bin/env nemesis
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

"""Script to run pylith Python test suite.

Run `coverage report` to generate a report (included).
Run `coverage html -d DIR` to generate an HTML report in directory `DIR`.
"""


import unittest
import sys
import importlib

DIRS = [
    "apps",
    "bc",
    # "faults",
    # "friction",
    # "materials",
    # "meshio",
    # "mpi",
    # "problems",
    # "tests",
    # "topology",
    # "utils",
]


class TestApp(object):
    """Application to run tests.
    """
    cov = None
    try:
        import coverage
        src_dirs = [
            "pylith.apps",
            "pylith.bc",
            "pylith.faults",
            "pylith.friction",
            "pylith.materials",
            "pylith.meshio",
            "pylith.mpi",
            "pylith.problems",
            "pylith.tests",
            "pylith.topology",
            "pylith.utils",
        ]
        cov = coverage.Coverage(source=src_dirs)
    except ImportError:
        pass

    def main(self):
        """
        Run the application.
        """
        if self.cov:
            self.cov.start()

        success = unittest.TextTestRunner(
            verbosity=2).run(self._suite()).wasSuccessful()
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
        for d in DIRS:
            sys.path.append(f"./{d}")

        test_cases = []
        for d in DIRS:
            mod = importlib.import_module(d)
            test_cases += mod.test_classes()

        suite = unittest.TestSuite()
        for test_case in test_cases:
            suite.addTest(unittest.makeSuite(test_case))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestApp().main()


# End of file
