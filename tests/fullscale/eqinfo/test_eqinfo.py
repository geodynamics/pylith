#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
from pylith.testing.FullTestApp import TestDriver


import TestEqInfoLine
import TestEqInfoTri
import TestEqInfoQuad


TEST_CLASSES = (
    *TestEqInfoLine.test_classes(),
    *TestEqInfoTri.test_classes(),
    *TestEqInfoQuad.test_classes(),
)



# ----------------------------------------------------------------------
if __name__ == '__main__':
    TestDriver(test_classes=TEST_CLASSES).run()


# End of file
