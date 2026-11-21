# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


from . import (
    TestCollectVersionInfo,
    TestConstants,
    TestEmptyBin,
    TestNullComponent,
    TestDumpParameters,
    TestDumpParametersAscii,
    TestDumpParametersJson,
    TestDependenciesVersion,
    TestPylithVersion,
    TestProfiling,
)


def test_modules():
    modules = [
        TestCollectVersionInfo,
        TestConstants,
        TestEmptyBin,
        TestNullComponent,
        TestDumpParameters,
        TestDumpParametersAscii,
        TestDumpParametersJson,
        TestDependenciesVersion,
        TestPylithVersion,
        TestProfiling,
    ]
    return modules


# End of file
