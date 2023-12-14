# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_linearelasticity_nofaults2d

TESTS = \
	run_UniformStrain2D.sh \
	run_Gravity2D.sh \
	run_GravityRefState2D.sh \
	run_BodyForce2D.sh \
	run_RigidBodyAcc2D.sh \
	run_PlanePWave2D.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = $(MMS_DRIVER)

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

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

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


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
