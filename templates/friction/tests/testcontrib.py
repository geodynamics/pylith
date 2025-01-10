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

import unittest


def suite():

    suite = unittest.TestSuite()

    from TestViscousFriction import TestViscousFriction
    suite.addTest(unittest.makeSuite(TestViscousFriction))

    return suite


def main():
    unittest.TextTestRunner(verbosity=2).run(suite())
    return


if __name__ == '__main__':
    main()


# End of file
