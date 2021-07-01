#!/usr/bin/env nemesis
#
# ======================================================================
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
# ======================================================================
#

import unittest


def suite():

    suite = unittest.TestSuite()

    from TestPlaneStrainState import TestPlaneStrainState
    suite.addTest(unittest.makeSuite(TestPlaneStrainState))

    return suite


def main():
    unittest.TextTestRunner(verbosity=2).run(suite())
    return


if __name__ == '__main__':
    main()


# End of file
