# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = elasticity
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = test_isotropiclinearelasticity

check_PROGRAMS = $(TESTS)


test_isotropiclinearelasticity_SOURCES = \
	TestElasticity.cc \
	TestIsotropicLinearElasticity.cc \
	TestIsotropicLinearElasticity2D_UniformStrain.cc \
	test_material.cc


noinst_HEADERS = 

dist_noinst_DATA = \
	data/tri_onecell.mesh \
	data/tri_fourcells.mesh \
	data/tri.mesh \
	data/quad_onecell.mesh \
	data/quad_fourcells.mesh \
	data/quad_small.mesh \
	data/quad_aligned.mesh \
	data/tet_small.mesh \
	data/hex.mesh


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


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/mmstests/elasticty/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
