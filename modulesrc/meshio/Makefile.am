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

subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = meshiomodule.la

meshiomodule_la_LDFLAGS = -module

meshiomodule_la_SOURCES = meshio.pyxe

nodist_meshiomodule_la_SOURCES = \
	meshio.c meshio_embed.cpp meshio_embed.h

meshiomodule_la_LIBADD = \
	$(top_builddir)/libsrc/meshio/libpylithmeshio.la

INCLUDES += -I$(PYTHON_INCDIR)

meshio.pyx meshio_embed.cpp  meshio_embed.h: meshio.pyxe
	cp $< . && pyrexembed.py meshio.pyxe && rm -f meshio.pyxe
meshio_embed.cpp: meshio_embed.h
meshio_embed.h: meshio.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = meshio.pyx meshio.c *_embed.*

# End of file 
