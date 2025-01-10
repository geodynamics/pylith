# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
	developer \
	docs \
	release-notes \
	templates \
	share

dist_noinst_PYTHON = \
	applications/pylith \
	applications/pylith_cfgsearch \
	applications/pylith_dumpparameters \
	applications/pylith_eqinfo \
	applications/pylith_genxdmf \
	applications/pylith_runner \
	applications/pylith_viz \
	applications/pylith_powerlaw_gendb


EXTRA_DIST = \
	AUTHORS \
	CITATION.cff \
	CHANGES.md \
	DEPENDENCIES \
	INSTALL \
	LICENSE.md \
	NEWS \
	README.md \
	setup.cfg \
	setup.py \
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
