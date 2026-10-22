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

from pylith.testing.FullTestApp import TestDriver, FullTestCase

import TestAxialStrainRateGenMaxwell
import TestAxialTractionMaxwell
import TestAxialStrainGenMaxwell

TEST_CLASSES = [
    *TestAxialStrainRateGenMaxwell.test_classes(),
    *TestAxialTractionMaxwell.test_classes(),
    *TestAxialStrainGenMaxwell.test_classes(),
]

# ----------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    TestDriver(test_classes=TEST_CLASSES).run()


# End of file
