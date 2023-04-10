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

include $(top_srcdir)/tests/check_cppunit.am

TESTS = \
	test_uniformstrain3d \
	test_gravity3d \
	test_gravityrefstate3d \
	test_bodyforce3d

check_PROGRAMS = $(TESTS)


test_uniformstrain3d_SOURCES = \
		TestLinearElasticity.cc \
		UniformStrain3D.cc \
		$(top_srcdir)/tests/src/driver_cppunit.cc

test_gravity3d_SOURCES = \
		TestLinearElasticity.cc \
		Gravity3D.cc \
		$(top_srcdir)/tests/src/driver_cppunit.cc

test_gravityrefstate3d_SOURCES = \
		TestLinearElasticity.cc \
		GravityRefState3D.cc \
		$(top_srcdir)/tests/src/driver_cppunit.cc

test_bodyforce3d_SOURCES = \
		TestLinearElasticity.cc \
		BodyForce3D.cc \
		$(top_srcdir)/tests/src/driver_cppunit.cc


dist_noinst_HEADERS = \
	TestLinearElasticity.hh	


dist_noinst_DATA = \
	data/tet.mesh \
	data/tet.msh \
	data/hex.mesh \
	data/hex.msh

noinst_TMP =

export_datadir = $(abs_builddir)/data


# End of file
