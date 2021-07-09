# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

subpackage = elasticity
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

TESTS = \
	test_isotropiclinearelasticity2d \
	test_isotropiclinearelasticity3d_uniformstrain \
	test_isotropiclinearelasticity3d_gravity \
	test_isotropiclinearelasticity3d_gravityrefstate \
	test_isotropiclinearelasticity3d_bodyforce

check_PROGRAMS = $(TESTS)


test_isotropiclinearelasticity2d_SOURCES = \
		TestElasticity.cc \
		TestIsotropicLinearElasticity.cc \
		TestIsotropicLinearElasticity2D_UniformStrain.cc \
		TestIsotropicLinearElasticity2D_Gravity.cc \
		TestIsotropicLinearElasticity2D_GravityRefState.cc \
		TestIsotropicLinearElasticity2D_BodyForce.cc \
		test_driver.cc

test_isotropiclinearelasticity3d_uniformstrain_SOURCES = \
		TestElasticity.cc \
		TestIsotropicLinearElasticity.cc \
		TestIsotropicLinearElasticity3D_UniformStrain.cc \
		test_driver.cc

test_isotropiclinearelasticity3d_gravity_SOURCES = \
		TestElasticity.cc \
		TestIsotropicLinearElasticity.cc \
		TestIsotropicLinearElasticity3D_Gravity.cc \
		test_driver.cc

test_isotropiclinearelasticity3d_gravityrefstate_SOURCES = \
		TestElasticity.cc \
		TestIsotropicLinearElasticity.cc \
		TestIsotropicLinearElasticity3D_GravityRefState.cc \
		test_driver.cc

test_isotropiclinearelasticity3d_bodyforce_SOURCES = \
		TestElasticity.cc \
		TestIsotropicLinearElasticity.cc \
		TestIsotropicLinearElasticity3D_BodyForce.cc \
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
	data/quad_distorted.mesh \
	data/tet.mesh \
	data/hex.mesh \
	data/tri_twocells.mesh \
	data/tri_fourcells.mesh \
	data/quad_fourcells.mesh

noinst_TMP =

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/mmstests/elasticity/data
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then $(MKDIR_P) $(export_datadir);  for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file
