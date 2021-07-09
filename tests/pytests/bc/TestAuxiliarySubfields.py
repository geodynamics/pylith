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
# @file tests/pytests/bc/TestAuxiliarySubfields.py
#
# @brief Unit testing of Python AuxSubfieldsTimeDependent and AuxSubfieldsAbsorbingDampers object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
import pylith.bc.AuxSubfieldsTimeDependent
import pylith.bc.AuxSubfieldsAbsorbingDampers


class TestAuxSubfieldsTimeDependent(TestComponent):
    """Unit testing of AuxSubfieldsTimeDependent object.
    """
    _class = pylith.bc.AuxSubfieldsTimeDependent.AuxSubfieldsTimeDependent
    _factory = pylith.bc.AuxSubfieldsTimeDependent.auxiliary_subfields


class TestAuxSubfieldsAbsorbingDampers(TestComponent):
    """Unit testing of AuxSubfieldsAbsorbingDampers object.
    """
    _class = pylith.bc.AuxSubfieldsAbsorbingDampers.AuxSubfieldsAbsorbingDampers
    _factory = pylith.bc.AuxSubfieldsAbsorbingDampers.auxiliary_subfields


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [TestAuxSubfieldsTimeDependent, TestAuxSubfieldsAbsorbingDampers]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
