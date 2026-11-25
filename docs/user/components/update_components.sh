# -*- sh -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


~/src/cig/pythia/bin/pyre_doc_components.py --package=PyLith --src-path=$PYLITH_DIR/lib/python$PYTHON_VERSION/site-packages/pylith --output-path=./ --manual-toc=toc.json
cat index.md | sed -e "s#apps/index.md#implementations.md\\napps/index.md#g" > tmp.md && mv -f tmp.md index.md

# Must edit PyLithApp.md and ConvertMeshApp.md manually.


# End of file
