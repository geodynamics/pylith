#!/bin/bash
# -*- sh -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


pylith elastic.cfg elastic_dt01.cfg
pylith elastic.cfg elastic_dt02.cfg
pylith elastic.cfg elastic_dt05.cfg
pylith elastic.cfg elastic_dt10.cfg

pylith plastic.cfg plastic_dt01.cfg
pylith plastic.cfg plastic_dt02.cfg
pylith plastic.cfg plastic_dt05.cfg
pylith plastic.cfg plastic_dt10.cfg

./plot_invars.py


# End of file
