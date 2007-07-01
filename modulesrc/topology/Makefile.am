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

subpkgpyexec_LTLIBRARIES = topologymodule.la

topologymodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

topologymodule_la_SOURCES = topology.pyxe

nodist_topologymodule_la_SOURCES = \
	topology.c topology_embed.cpp topology_embed.h

topologymodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)

if ENABLE_CUBIT
  topologymodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
  topologymodule_la_LIBADD += $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

PYTHON_EGG_PYXFLAGS = $(PYTHON_EGG_CPPFLAGS)/mpi

topology.pyx topology_embed.cpp  topology_embed.h: topology.pyxe
	pyrexembed topology.pyxe
topology.pyxe: $(srcdir)/topology.pyxe.src
	cp $(srcdir)/topology.pyxe.src $@
topology_embed.cpp: topology_embed.h
topology_embed.h: topology.pyx

.pyx.c:
	pyrexc $(PYTHON_EGG_PYXFLAGS) $<

CLEANFILES = topology.pyxe topology.pyx topology.c *_embed.*

# End of file 
