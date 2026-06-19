# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

subpackage = utils
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = \
	_utils.la

subpkgpython_PYTHON = \
	utils.py

utils_swig_sources = \
	utils.i \
	../include/chararray.i \
	pylith_general.i \
	PyreComponent.i \
	PylithVersion.i \
	DependenciesVersion.i \
	TestArray.i \
	constants.i \
	PetscVersion.i \
	Application.i \
	PetscOptions.i \
	EventLogger.i

utils_swig_generated = \
	utils_wrap.cxx \
	utils.py

swig_sources = $(utils_swig_sources)

swig_generated = $(utils_swig_generated)

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
$(srcdir)/utils_wrap.cxx $(srcdir)/utils.py: $(utils_swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/utils_wrap.cxx $(srcdir)/utils.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/utils_wrap.cxx \
	$(srcdir)/utils.py


# End of file
