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

include $(top_srcdir)/tests/check_catch2.am


TESTS = \
	run_uniformstrain3d.sh \
	run_gravity3d.sh \
	run_gravityrefstate3d.sh \
	run_bodyforce3d.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = mmstest_linearelasticity_nofaults3d

run_%3d.sh: mmstest_linearelasticity_nofaults3d

mmstest_linearelasticity_nofaults3d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestLinearElasticity.cc \
	TestCases.cc \
	UniformStrain3D.cc \
	Gravity3D.cc \
	GravityRefState3D.cc \
	BodyForce3D.cc 


dist_noinst_HEADERS = \
	TestLinearElasticity.hh	 \
	UniformStrain3D.hh \
	Gravity3D.hh \
	GravityRefState3D.hh \
	BodyForce3D.hh


dist_noinst_DATA = \
	data/tet.mesh \
	data/tet.msh \
	data/hex.mesh \
	data/hex.msh

noinst_TMP =

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
