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

TESTS = mmstest_linearelasticity_nofaults2d 

check_PROGRAMS = $(TESTS)


mmstest_linearelasticity_nofaults2d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestLinearElasticity.cc \
	TestCases.cc \
	UniformStrain2D.cc \
	Gravity2D.cc \
	GravityRefState2D.cc \
	BodyForce2D.cc \
	RigidBodyAcc2D.cc \
	PlanePWave2D.cc 



dist_noinst_HEADERS = \
	TestLinearElasticity.hh	\
	UniformStrain2D.hh \
	Gravity2D.hh \
	GravityRefState2D.hh \
	BodyForce2D.hh \
	RigidBodyAcc2D.hh \
	PlanePWave2D.hh


dist_noinst_DATA = \
	data/tri.mesh \
	data/tri.msh \
	data/tri_twocells.mesh \
	data/tri_fourcells.mesh \
	data/quad.mesh \
	data/quad.msh \
	data/quad_distorted.mesh \
	data/quad_fourcells.mesh


# 'export' the input files by performing a mock install
export_datadir = $(abs_builddir)/data


# End of file
