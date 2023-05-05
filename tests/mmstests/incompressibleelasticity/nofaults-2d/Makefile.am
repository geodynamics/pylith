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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_incompressibleelasticity_nofaults2d

TESTS = \
	run_UniformPressure2D.sh \
	run_UniformShear2D.sh \
	run_Gravity2D.sh \
	run_BodyForce2D.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = $(MMS_DRIVER)


mmstest_incompressibleelasticity_nofaults2d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestIncompressibleElasticity.cc \
	TestCases.cc \
	UniformPressure2D.cc \
	UniformShear2D.cc \
	Gravity2D.cc \
	BodyForce2D.cc 

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
	TestIncompressibleElasticity.hh \
	UniformPressure2D.hh \
	UniformShear2D.hh \
	Gravity2D.hh \
	BodyForce2D.hh

dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh \
	data/tet.mesh \
	data/hex.mesh

noinst_TMP =

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
