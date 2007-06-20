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

meshiomodule_la_LDFLAGS = -module -no-undefined -avoid-version \
	$(PYTHON_LA_LDFLAGS)

meshiomodule_la_SOURCES = meshio.pyxe

nodist_meshiomodule_la_SOURCES = \
	meshio.c meshio_embed.cpp meshio_embed.h

meshiomodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  meshiomodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

PYREX_INCLUDES = -I$(top_srcdir)/modulesrc/topology

meshio.pyx meshio_embed.cpp  meshio_embed.h: meshio.pyxe
	pyrexembed meshio.pyxe
meshio.pyxe: $(srcdir)/meshio.pyxe.src
	cp $(srcdir)/meshio.pyxe.src $@
meshio_embed.cpp: meshio_embed.h
meshio_embed.h: meshio.pyx

.pyx.c:
	pyrexc $< $(PYREX_INCLUDES)

CLEANFILES = meshio.pyxe meshio.pyx meshio.c *_embed.*

# End of file 
