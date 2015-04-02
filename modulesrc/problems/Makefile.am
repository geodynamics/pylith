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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = problems
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = _problemsmodule.la

subpkgpyexec_PYTHON = problems.py

swig_sources = \
	problems.i \
	Formulation.i \
	Explicit.i \
	Implicit.i \
	Solver.i \
	SolverLinear.i \
	SolverNonlinear.i \
	SolverLumped.i


swig_generated = \
	problems_wrap.cxx \
	problems.py

_problemsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__problemsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_problemsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _problemsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_problemsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/problems_wrap.cxx $(srcdir)/problems.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/problems_wrap.cxx \
	$(srcdir)/problems.py


# End of file 
