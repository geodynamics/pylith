# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = friction
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _friction.la

subpkgpython_PYTHON = friction.py

swig_sources = \
	friction.i \
	../utils/PyreComponent.i

swigincludedir = $(pkgdatadir)/swig/$(subpackage)
#swiginclude_HEADERS = \
#	FrictionModel.i

swig_generated = \
	friction_wrap.cxx \
	friction.py

_friction_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__friction_la_SOURCES = $(swig_sources) $(swig_generated)

_friction_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _friction_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/friction_wrap.cxx $(srcdir)/friction.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/friction_wrap.cxx $(srcdir)/friction.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/friction_wrap.cxx \
	$(srcdir)/friction.py


# End of file 
