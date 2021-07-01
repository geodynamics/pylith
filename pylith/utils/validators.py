# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
"""Python module for some additional Pyre validation.
"""

def notEmptyList(value):
    print(f"LENGTH of {value}= {len(value)}")
    if len(value) > 0:
        return value
    raise ValueError("Nonempty list required.")

def notEmptyString(value):
    if len(value) > 0:
        return value
    raise ValueError("Nonempty string required.")


# End of file
