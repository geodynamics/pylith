# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = problems
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _problems.la

subpkgpython_PYTHON = problems.py

swig_sources = \
	problems.i \
	../utils/PyreComponent.i \
	../topology/FieldBase.i \
	../include/physicsarray.i \
	../include/outputarray.i \
	../include/scalartypemaps.i \
	Problem.i \
	TimeDependent.i \
	GreensFns.i \
	Physics.i \
	ObserverSoln.i \
	ObserverPhysics.i \
	InitialCondition.i \
	InitialConditionDomain.i \
	InitialConditionPatch.i \
	ProgressMonitor.i \
	ProgressMonitorTime.i


swig_generated = \
	problems_wrap.cxx \
	problems.py

_problems_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__problems_la_SOURCES = $(swig_sources) $(swig_generated)

_problems_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _problems_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/problems_wrap.cxx $(srcdir)/problems.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/problems_wrap.cxx $(srcdir)/problems.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/problems_wrap.cxx \
	$(srcdir)/problems.py


# End of file
