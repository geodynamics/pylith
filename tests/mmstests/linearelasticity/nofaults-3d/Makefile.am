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
	test_uniformstrain3d \
	test_gravity3d \
	test_gravityrefstate3d \
	test_bodyforce3d

check_PROGRAMS = $(TESTS)


test_uniformstrain3d_SOURCES = \
		TestLinearElasticity.cc \
		UniformStrain3D.cc \
		test_driver.cc

test_gravity3d_SOURCES = \
		TestLinearElasticity.cc \
		Gravity3D.cc \
		test_driver.cc

test_gravityrefstate3d_SOURCES = \
		TestLinearElasticity.cc \
		GravityRefState3D.cc \
		test_driver.cc

test_bodyforce3d_SOURCES = \
		TestLinearElasticity.cc \
		BodyForce3D.cc \
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
	data/tet.mesh \
	data/tet.msh \
	data/hex.mesh \
	data/hex.msh

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
