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
# Copyright (c) 2010-2015 University of California, Davis
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
	DEPENDENCIES

# End of file 
