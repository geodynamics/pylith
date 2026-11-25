# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


__all__ = [
    "FullTestApp",
    "UnitTestApp",
]


# ----------------------------------------------------------------------
def has_h5py():
    if not "flag" in dir(has_h5py):
        try:
            import h5py

            has_h5py.flag = True
        except ImportError:
            print("WARNING: Cannot find h5py Python modele.")
            print("         Tests limited to running PyLith without errors.")
            print("         Install h5py (available via the installer utility) ")
            print("         in order to enable verification of output.")
            has_h5py.flag = False
    return has_h5py.flag


# End of file
