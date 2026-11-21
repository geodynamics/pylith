#!/bin/bash
# -*- sh -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


for fault in main east west; do
    cat faulttrace_${fault}_latlon.txt | cs2cs epsg:4326 epsg:32611 > faulttrace_${fault}_utm.txt
done


# End of file
