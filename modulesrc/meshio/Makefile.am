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
	$(top_builddir)/libsrc/meshio/libpylithmeshio.la $(PETSC_LIB)

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

meshio.pyx meshio_embed.cpp  meshio_embed.h: meshio.pyxe
	if [ "${VPATH}" != "" ]; then cp $< .; fi && pyrexembed meshio.pyxe && if [ "${VPATH}" != "" ]; then rm -f meshio.pyxe; fi
meshio_embed.cpp: meshio_embed.h
meshio_embed.h: meshio.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = meshio.pyx meshio.c *_embed.*

# End of file 
