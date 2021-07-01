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

subpackage = feassemble
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _feassemble.la

subpkgpython_PYTHON = feassemble.py

swig_sources = \
	feassemble.i

swig_generated = \
	feassemble_wrap.cxx \
	feassemble.py

_feassemble_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__feassemble_la_SOURCES = $(swig_sources) $(swig_generated)

_feassemble_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _feassemble_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/feassemble_wrap.cxx $(srcdir)/feassemble.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/feassemble_wrap.cxx $(srcdir)/feassemble.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/feassemble_wrap.cxx \
	$(srcdir)/feassemble.py


# End of file
