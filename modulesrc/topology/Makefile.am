# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
	Fields.i \
	SolutionFields.i

swig_generated = \
	topology_wrap.cxx \
	topology.py

_topologymodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__topologymodule_la_SOURCES = $(swig_sources) $(swig_generated)

_topologymodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if NO_UNDEFINED
_topologymodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/topology_wrap.cxx $(srcdir)/topology.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/topology_wrap.cxx \
	$(srcdir)/topology.py


# End of file 
