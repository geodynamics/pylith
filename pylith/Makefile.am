# -*- Makefile -*-
#
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
#

if ENABLE_CUBIT
make-manifest:
else
make-manifest:
	echo "exclude meshio/MeshIOCubit.py" > $@
endif

install-exec-local: make-manifest
	$(PYTHON) -m pip install $(top_srcdir)

# End of file
