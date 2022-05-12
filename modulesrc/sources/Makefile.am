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

subpackage = sources
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _sources.la

subpkgpython_PYTHON = sources.py

swig_sources = \
	sources.i \
	../problems/Physics.i \
	Source.i \
	PointForce.i \
	SourceTimeFunctionPointForce.i \
	RickerWavelet.i \
	GaussianWavelet.i \	
	WellboreSource.i \
	../utils/PyreComponent.i

swigincludedir = $(pkgdatadir)/swig/$(subpackage)
swiginclude_HEADERS = \
	Source.i

swig_generated = \
	sources_wrap.cxx \
	sources.py

_sources_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__sources_la_SOURCES = $(swig_sources) $(swig_generated)

_sources_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _sources_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/sources_wrap.cxx $(srcdir)/sources.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/sources_wrap.cxx $(srcdir)/sources.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/sources_wrap.cxx \
	$(srcdir)/sources.py


# End of file