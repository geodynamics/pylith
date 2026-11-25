# -*- makefile -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


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
