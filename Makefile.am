# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = \
	libsrc \
	modulesrc \
	pylith \
	applications \
	templates \
	share

if ENABLE_TESTING
  SUBDIRS += \
	unittests 
endif
if ENABLE_FULL_TESTING
  SUBDIRS += \
	tests_auto
endif

DIST_SUBDIRS = $(SUBDIRS) \
	examples \
	tests \
	doc

if ENABLE_DOCUMENTATION
  SUBDIRS += doc
endif

EXTRA_DIST = \
	CHANGES \
	DEPENDENCIES \
	setup.py \
	archimedes/__init__.py

install-data-hook:
	top_builddir=`$(am__cd) $(top_builddir) && pwd` && \
		cd $(top_srcdir) && \
		$(PYTHON) setup.py install_deps -H None -f $$top_builddir/deps --prefix=$(prefix) && \
		$(PYTHON) setup.py egg_info -e $(pythondir)


# End of file 
