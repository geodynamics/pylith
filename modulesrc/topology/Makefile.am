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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = topology
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _topologymodule.la

subpkgpyexec_PYTHON = topology.py

swig_sources = \
	topology.i \
	Mesh.i \
	SubMesh.i \
	MeshOps.i \
	FieldBase.i \
	Field.i \
	MultiField.i \
	Fields.i \
	PackedFields.i \
	SolutionFields.i \
	Jacobian.i \
	Distributor.i \
	RefineUniform.i \
	ReverseCuthillMcKee.i

swig_generated = \
	topology_wrap.cxx \
	topology.py

_topologymodule_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__topologymodule_la_SOURCES = $(swig_sources) $(swig_generated)

_topologymodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if NO_UNDEFINED
_topologymodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += -g $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_CC_INCLUDES)

if ENABLE_SWIG
$(srcdir)/topology_wrap.cxx $(srcdir)/topology.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/topology_wrap.cxx \
	$(srcdir)/topology.py


# End of file 
