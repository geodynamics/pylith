# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Python module for some additional Pyre validation.
"""

def notEmptyList(value):
    if len(value) > 0:
        return value
    raise ValueError("Nonempty list required.")

def notEmptyString(value):
    if len(value) > 0:
        return value
    raise ValueError("Nonempty string required.")


# End of file
