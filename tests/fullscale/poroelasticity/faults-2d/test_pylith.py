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

import TestPressureGradient

TEST_CLASSES = [
    *TestPressureGradient.test_classes(),
]

# ----------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    TestDriver(test_classes=TEST_CLASSES).run()


# End of file
