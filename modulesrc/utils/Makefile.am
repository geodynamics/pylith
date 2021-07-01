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

subpackage = utils
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = \
	_petsc.la \
	_utils.la

subpkgpython_PYTHON = \
	petsc.py \
	utils.py

petsc_swig_sources = \
	petsc.i \
	petsc_general.i \
	petsc_memory.i \
	petsc_mat.i

petsc_swig_generated = \
	petsc_wrap.cxx \
	petsc.py

utils_swig_sources = \
	utils.i \
	pylith_general.i \
	PylithVersion.i \
	PetscVersion.i \
	DependenciesVersion.i \
	EventLogger.i \
	PyreComponent.i \
	TestArray.i \
	constdefs.i

utils_swig_generated = \
	utils_wrap.cxx \
	utils.py

swig_sources = $(petsc_swig_sources) $(utils_swig_sources)

swig_generated = $(petsc_swig_generated) $(utils_swig_generated)

# PETSc module
_petsc_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__petsc_la_SOURCES = $(petsc_swig_sources) $(petsc_swig_generated)

_petsc_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

# utils module
_utils_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__utils_la_SOURCES = $(utils_swig_sources) $(utils_swig_generated)

_utils_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/petsc_wrap.cxx $(srcdir)/petsc.py: $(petsc_swig_sources)
	$(SWIG) -Wall -c++ -python $<

$(srcdir)/utils_wrap.cxx $(srcdir)/utils.py: $(utils_swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/petsc_wrap.cxx $(srcdir)/petsc.py $(srcdir)/utils_wrap.cxx $(srcdir)/utils.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/petsc_wrap.cxx \
	$(srcdir)/petsc.py \
	$(srcdir)/utils_wrap.cxx \
	$(srcdir)/utils.py


# End of file
