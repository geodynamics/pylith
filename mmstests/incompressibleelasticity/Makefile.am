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

subpackage = incompressibleelasticity
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = test_isotropiclinearelasticity

check_PROGRAMS = $(TESTS)


test_isotropiclinearelasticity_SOURCES = \
	TestIncompressibleElasticity.cc \
	TestIsotropicLinearIncompElasticity.cc \
	TestIsotropicLinearIncompElasticity2D_UniformPressure.cc \
	TestIsotropicLinearIncompElasticity2D_UniformShear.cc \
	TestIsotropicLinearIncompElasticity2D_Gravity.cc \
	test_driver.cc




noinst_HEADERS =


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
	data/quad_aligned.mesh \
	data/tet.mesh \
	data/hex.mesh

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/mmstests/incompressibleelasticity/data
export-data: $(dist_noinst_DATA)
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then mkdir $(export_datadir);  for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


leakcheck-isotropiclinearelasticity:
	valgrind --log-file=valgrind_isotropiclinearelasticity.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/test_isotropiclinearelasticity



# End of file
