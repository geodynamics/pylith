#!/usr/bin/env nemesis
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from pylith.testing.FullTestApp import TestDriver, FullTestCase

import TestAxialTractionMaxwell
import TestAxialStrainGenMaxwell
import TestAxialStrainRateGenMaxwell

TEST_CLASSES = [
    *TestAxialTractionMaxwell.test_classes(),
    *TestAxialStrainGenMaxwell.test_classes(),
    *TestAxialStrainRateGenMaxwell.test_classes(),
]

# ----------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    TestDriver(test_classes=TEST_CLASSES).run()


# End of file
