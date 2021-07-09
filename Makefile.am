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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = \
	libsrc \
	modulesrc \
	pylith

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
	DEPENDENCIES \
	docker/pylith-testenv \
	docker/pylith-docenv


if ENABLE_TEST_COVERAGE
coverage-%.info: 
	lcov --directory libsrc --capture --output-file $@
	lcov --remove $@ '/usr/*' '/Library/*' "v1" "*/testing/TestDriver.*" -o $@
	lcov --list $@

coverage-pytests.xml:
	$(PYTHON_COVERAGE) xml -o coverage-pytests.xml tests/pytests/coverage.xml

coverage-%: coverage-%.info
	genhtml -o $@ $^

coverage-pytests: 
	$(PYTHON_COVERAGE) html -d coverage-pytests tests/pytests/coverage.xml

clean-coverage:
	$(RM) `find . -name "*.gcda"`
	$(RM) `find . -name ".coverage"`

else
coverage-libtests.info coverage-mmstests.info coverage-pytests.xml coverage-fullscale.info coverage-libtests coverage-mmstests coverage-pytests coverage-fullscale:
	$(warning Test coverage not enabled. Ignoring test coverage targets)

clean-coverage:
	$(warning Test coverage not enabled. Nothing to clean.)

endif


# Remove .info files from list of intermediate files
.SECONDARY: coverage-libtests.info coverage-mmstests.info coverage-pytests.xml coverage-fullscale.info

# End of file
