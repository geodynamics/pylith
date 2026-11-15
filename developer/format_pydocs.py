#!/usr/bin/env python3
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import pathlib
import re

PATTERN_DEF = r"(def \w+.*:\n\s+\"\"\")\s?\n\s+(\w+)"
PATTERN_CLASS = r"(class \w+.*:\n\s+\"\"\")\s?\n\s+(\w+)"


def replacement(match):
    return "".join(match.groups())


files = pathlib.Path(".").rglob("*.py")
for filename in files:
    with open(filename, "r") as fin:
        lines = fin.read()
        newlines = re.sub(PATTERN_DEF, replacement, lines)
        newlines = re.sub(PATTERN_CLASS, replacement, newlines)
    with open(filename, "w") as fout:
        fout.write(newlines)


# End of file
