# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
