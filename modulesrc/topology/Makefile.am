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

topologymodule_la_LDFLAGS = -module

topologymodule_la_SOURCES = topology.pyxe

nodist_topologymodule_la_SOURCES = \
	topology.c topology_embed.cpp topology_embed.h

topologymodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB)

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

topology.pyx topology_embed.cpp  topology_embed.h: topology.pyxe
	pyrexembed topology.pyxe
topology.pyxe: $(srcdir)/topology.pyxe.src
	cp $(srcdir)/topology.pyxe.src $@
topology_embed.cpp: topology_embed.h
topology_embed.h: topology.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = topology.pyxe topology.pyx topology.c *_embed.*

# End of file 
