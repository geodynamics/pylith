# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

subpackage = petsc
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = \
	_petsc.la

subpkgpython_PYTHON = \
	petsc.py

petsc_swig_sources = \
	petsc.i \
	../include/chararray.i \
	PetscVersion.i \
	Application.i \
	Options.i \
	EventLogger.i

petsc_swig_generated = \
	petsc_wrap.cxx \
	petsc.py

swig_sources = $(petsc_swig_sources)

swig_generated = $(petsc_swig_generated)

# petsc module
_petsc_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__petsc_la_SOURCES = $(petsc_swig_sources) $(petsc_swig_generated)

_petsc_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/petsc_wrap.cxx $(srcdir)/petsc.py: $(petsc_swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/petsc_wrap.cxx $(srcdir)/petsc.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/petsc_wrap.cxx \
	$(srcdir)/petsc.py


# End of file
