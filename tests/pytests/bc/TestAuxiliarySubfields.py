#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
