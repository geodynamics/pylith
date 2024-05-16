# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_linearelasticity_nofaults3d

TESTS = \
	run_UniformStrain3D.sh \
	run_Gravity3D.sh \
	run_GravityRefState3D.sh \
	run_BodyForce3D.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearelasticity_nofaults3d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestLinearElasticity.cc \
	TestCases.cc \
	UniformStrain3D.cc \
	Gravity3D.cc \
	GravityRefState3D.cc \
	BodyForce3D.cc 

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
	TestLinearElasticity.hh	 \
	UniformStrain3D.hh \
	Gravity3D.hh \
	GravityRefState3D.hh \
	BodyForce3D.hh


dist_noinst_DATA = \
	data/tet.msh \
	data/hex.msh

noinst_TMP =

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
