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

subpackage = poroelasticity
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = \
	test_isotropiclinearporoelasticity2d \
	test_isotropiclinearporoelasticity2d_QS_LT \
	test_isotropiclinearporoelasticity2d_QS_TT \
	test_isotropiclinearporoelasticity3d \
	test_isotropiclinearporoelasticity3d_QS_LT \
	test_isotropiclinearporoelasticity3d_QS_TT 	

check_PROGRAMS = $(TESTS)

test_isotropiclinearporoelasticity2d_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity2D_QS_LT.cc \
		TestIsotropicLinearPoroelasticity2D_QS_TT.cc \
		test_driver.cc

test_isotropiclinearporoelasticity2d_QS_LT_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity2D_QS_LT.cc \
		test_driver.cc

test_isotropiclinearporoelasticity2d_QS_TT_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity2D_QS_TT.cc \
		test_driver.cc

test_isotropiclinearporoelasticity3d_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity3D_QS_LT.cc \
		TestIsotropicLinearPoroelasticity3D_QS_TT.cc \
		test_driver.cc

test_isotropiclinearporoelasticity3d_QS_LT_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity3D_QS_LT.cc \
		test_driver.cc

test_isotropiclinearporoelasticity3d_QS_TT_SOURCES = \
		TestPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity.cc \
		TestIsotropicLinearPoroelasticity3D_QS_TT.cc \
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
	data/quad.mesh \
	data/tet.mesh \
	data/hex.mesh 

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/mmstests/poroelasticity/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then $(MKDIR_P) $(export_datadir);  for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
