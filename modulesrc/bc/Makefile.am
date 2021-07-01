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

subpackage = bc
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _bc.la

subpkgpython_PYTHON = bc.py

swig_sources = \
	bc.i \
	../utils/PyreComponent.i \
	../problems/Physics.i \
	BoundaryCondition.i \
	DirichletTimeDependent.i \
	NeumannTimeDependent.i \
	AbsorbingDampers.i

swig_generated = \
	bc_wrap.cxx \
	bc.py

_bc_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__bc_la_SOURCES = $(swig_sources) $(swig_generated)

_bc_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _bc_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/bc_wrap.cxx $(srcdir)/bc.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/bc_wrap.cxx $(srcdir)/bc.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/bc_wrap.cxx \
	$(srcdir)/bc.py


# End of file
