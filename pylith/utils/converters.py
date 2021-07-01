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
"""Python module for simple conversions.
"""

import re


def string_to_list(option):
    """Convert string to list.

    List may be enclosed in [], {}, (), or not enclosed. List entries must be comma delimited.

    Args:
        option (str)
            String to convert.
    Returns: (list)
        List of strings.
    """
    PATTERN = r"[\{\[\(](.*)[\}\]\)]"

    match = re.search(PATTERN, option, re.DOTALL)
    optionStr = match.groups()[0] if match else option
    return [name.strip() for name in optionStr.split(",")]


# End of file
