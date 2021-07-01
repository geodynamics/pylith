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
		cd $(top_srcdir) && \
			$(PYTHON) setup.py egg_info && \
			PYTHONUSERBASE=$(DESTDIR)$(prefix) $(PYTHON) setup.py install --record record.txt --single-version-externally-managed --user && \
			cd $(abs_builddir)

# Pip install is too slow, because it makes a copy of the entire source tree, including .git.
# See https://github.com/pypa/pip/issues/2195
#PYTHONUSERBASE=$(DESTDIR)$(prefix) $(PYTHON) -m pip install --no-index --user $(top_srcdir)

# End of file
