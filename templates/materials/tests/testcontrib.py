#!/usr/bin/env nemesis
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


import unittest


def suite():

    suite = unittest.TestSuite()

    from TestPlaneStrainState import TestPlaneStrainState

    suite.addTest(unittest.makeSuite(TestPlaneStrainState))

    return suite


def main():
    unittest.TextTestRunner(verbosity=2).run(suite())
    return


if __name__ == "__main__":
    main()


# End of file
