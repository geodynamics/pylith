# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
