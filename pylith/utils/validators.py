# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


"""Python module for some additional Pyre validation."""


def notEmptyList(value):
    if len(value) > 0:
        return value
    raise ValueError("Nonempty list required.")


def notEmptyString(value):
    if len(value) > 0:
        return value
    raise ValueError("Nonempty string required.")


# End of file
