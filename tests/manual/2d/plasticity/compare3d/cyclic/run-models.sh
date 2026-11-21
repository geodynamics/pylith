# -*- sh -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


pylith onecell_2d.cfg quasistatic.cfg cycle_cont2d.cfg
pylith onecell_3d.cfg quasistatic.cfg cycle_cont3d.cfg

pylith onecell_2d.cfg quasistatic.cfg sin2d.cfg sin2d_quasistatic.cfg
pylith onecell_2d.cfg dynamic.cfg sin2d.cfg sin2d_dynamic.cfg
pylith onecell_3d.cfg quasistatic.cfg sin3d.cfg sin3d_quasistatic.cfg
pylith onecell_3d.cfg dynamic.cfg sin3d.cfg sin3d_dynamic.cfg


# End of file
