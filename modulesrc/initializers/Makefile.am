# -*- makefile -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


subpackage = initializers
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _initializers.la
subpkgpython_PYTHON = initializers.py

swig_sources = \
	initializers.i \
	phasesarray.i \
	Initializer.i \
	MeshReader.i \
	MeshWriter.i \
	MeshReordering.i \
	MeshRefiner.i \
	MeshDistributor.i \
	MeshInsertInterfaces.i \
	InitializePhase.i

swig_generated = \
	initializers_wrap.cxx \
	initializers.py

_initializers_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__initializers_la_SOURCES = $(swig_sources) $(swig_generated)

_initializers_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/initializers_wrap.cxx $(srcdir)/initializers.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<
else
$(srcdir)/initializers_wrap.cxx $(srcdir)/initializers.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/initializers_wrap.cxx \
	$(srcdir)/initializers.py


# End of file
