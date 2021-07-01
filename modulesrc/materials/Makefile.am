# -*- Makefile -*-
#
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
#

subpackage = materials
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _materials.la

subpkgpython_PYTHON = materials.py

swig_sources = \
	materials.i \
	../problems/Physics.i \
	Material.i \
	Elasticity.i \
	RheologyElasticity.i \
	IsotropicLinearElasticity.i \
	IsotropicLinearMaxwell.i \
	IsotropicLinearGenMaxwell.i \
	IsotropicPowerLaw.i \
	IncompressibleElasticity.i \
	RheologyIncompressibleElasticity.i \
	IsotropicLinearIncompElasticity.i \
	IsotropicLinearPoroelasticity.i \
	RheologyPoroelasticity.i \
	Poroelasticity.i \
	../utils/PyreComponent.i

swigincludedir = $(pkgdatadir)/swig/$(subpackage)
swiginclude_HEADERS = \
	Material.i

swig_generated = \
	materials_wrap.cxx \
	materials.py

_materials_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__materials_la_SOURCES = $(swig_sources) $(swig_generated)

_materials_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _materials_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/materials_wrap.cxx $(srcdir)/materials.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/materials_wrap.cxx $(srcdir)/materials.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/materials_wrap.cxx \
	$(srcdir)/materials.py


# End of file
