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

include $(top_srcdir)/check.am

TESTS = \
	test_uniformstrain2d \
	test_gravity2d \
	test_gravityrefstate2d \
	test_bodyforce2d \
	test_rigidbodyacc2d \
	test_planewave2d

check_PROGRAMS = $(TESTS)


test_uniformstrain2d_SOURCES = \
		TestLinearElasticity.cc \
		UniformStrain2D.cc \
		test_driver.cc

test_gravity2d_SOURCES = \
		TestLinearElasticity.cc \
		Gravity2D.cc \
		test_driver.cc

test_gravityrefstate2d_SOURCES = \
		TestLinearElasticity.cc \
		GravityRefState2D.cc \
		test_driver.cc

test_bodyforce2d_SOURCES = \
		TestLinearElasticity.cc \
		BodyForce2D.cc \
		test_driver.cc

test_rigidbodyacc2d_SOURCES = \
		TestLinearElasticity.cc \
		RigidBodyAcc2D.cc \
		test_driver.cc

test_planewave2d_SOURCES = \
		TestLinearElasticity.cc \
		PlanePWave2D.cc \
		test_driver.cc


dist_noinst_HEADERS = \
	TestLinearElasticity.hh	


AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_CC_INCLUDES)

AM_LDFLAGS = $(PYTHON_LA_LDFLAGS)

LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  LDADD += -lnetcdf
endif


dist_noinst_DATA = \
	data/tri.mesh \
	data/tri.msh \
	data/tri_twocells.mesh \
	data/tri_fourcells.mesh \
	data/quad.mesh \
	data/quad.msh \
	data/quad_distorted.mesh \
	data/quad_fourcells.mesh

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_builddir)/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then $(MKDIR_P) $(export_datadir);  for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
