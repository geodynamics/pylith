#!/usr/bin/env python3

import pathlib
import re

PATTERN_DEF = r"(def \w+.*\n\s+\"\"\")\s?\n\s+(\w+)"
PATTERN_CLASS = r"(class \w+.*\n\s+\"\"\")\s?\n\s+(\w+)"


def replacement(match):
    return "".join(match.groups())


files = pathlib.Path(".").glob("**/*.py")
for filename in files:
    with open(filename, 'r') as fin:
        lines = fin.read()
        newlines = re.sub(PATTERN_DEF, replacement, lines)
        newlines = re.sub(PATTERN_CLASS, replacement, newlines)
    with open(filename, 'w') as fout:
        fout.write(newlines)
