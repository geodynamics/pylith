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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = utils
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = \
	_petscmodule.la \
	_utilsmodule.la

subpkgpyexec_PYTHON = \
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
	EventLogger.i \
	TestArray.i \
	constdefs.i

utils_swig_generated = \
	utils_wrap.cxx \
	utils.py

swig_sources = $(petsc_swig_sources) $(utils_swig_sources)

swig_generated = $(petsc_swig_generated) $(utils_swig_generated)

# PETSc module
_petscmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__petscmodule_la_SOURCES = $(petsc_swig_sources) $(petsc_swig_generated)

_petscmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if NO_UNDEFINED
_petscmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

# utils module
_utilsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__utilsmodule_la_SOURCES = $(utils_swig_sources) $(utils_swig_generated)

_utilsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if NO_UNDEFINED
_utilsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/petsc_wrap.cxx $(srcdir)/petsc.py: $(petsc_swig_sources)
	$(SWIG) -Wall -c++ -python $<

$(srcdir)/utils_wrap.cxx $(srcdir)/utils.py: $(utils_swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/petsc_wrap.cxx \
	$(srcdir)/petsc.py \
	$(srcdir)/utils_wrap.cxx \
	$(srcdir)/utils.py


# End of file 
