# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = scales
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _scales.la
subpkgpython_PYTHON = scales.py

swig_sources = \
	scales.i \
	ScalesObj.i \
	ElasticityScales.i

swig_generated = \
	scales_wrap.cxx \
	scales.py

_scales_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__scales_la_SOURCES = $(swig_sources) $(swig_generated)

_scales_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/scales_wrap.cxx $(srcdir)/scales.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<
else
$(srcdir)/scales_wrap.cxx $(srcdir)/scales.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/scales_wrap.cxx \
	$(srcdir)/scales.py


# End of file
