# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = topology
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _topology.la

subpkgpython_PYTHON = topology.py

swig_sources = \
	topology.i \
	Mesh.i \
	MeshOps.i \
	FieldBase.i \
	Field.i \
	Distributor.i \
	RefineUniform.i \
	ReverseCuthillMcKee.i

swig_generated = \
	topology_wrap.cxx \
	topology.py

_topology_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__topology_la_SOURCES = $(swig_sources) $(swig_generated)

_topology_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/topology_wrap.cxx $(srcdir)/topology.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/topology_wrap.cxx $(srcdir)/topology.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/topology_wrap.cxx \
	$(srcdir)/topology.py


# End of file 
