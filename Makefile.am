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
# Copyright (c) 2010-2019 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = \
	libsrc \
	modulesrc \
	pylith \
	applications

if ENABLE_TESTING
  SUBDIRS += \
	tests
endif

DIST_SUBDIRS = $(SUBDIRS) \
	examples \
	doc \
	templates \
	share

EXTRA_DIST = \
	CHANGES \
	DEPENDENCIES

coverage-%.info: 
	lcov --directory libsrc --capture --output-file $@
	lcov --remove $@ '/usr/*' '/Library/*' "v1" "*/testing/TestDriver.*" -o $@
	lcov --list $@

coverage-%: coverage-%.info
	genhtml -o $@ $<

clean-coverage:
	$(RM) `find . -name "*.gcda"`

# Remove .info files from list of intermediate files
.SECONDARY: coverage-libtests.info coverage-mmstests.info coverage-fullscale.info

# End of file
