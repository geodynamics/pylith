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

subpackage = materials
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = materialsmodule.la

materialsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

materialsmodule_la_SOURCES = materials.pyxe

nodist_materialsmodule_la_SOURCES = \
	materials.c materials_embed.cpp materials_embed.h

materialsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  materialsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

materials.pyx materials_embed.cpp  materials_embed.h: materials.pyxe
	pyrexembed materials.pyxe
materials.pyxe: $(srcdir)/materials.pyxe.src
	cp $(srcdir)/materials.pyxe.src $@
materials_embed.cpp: materials_embed.h
materials_embed.h: materials.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = materials.pyxe materials.pyx materials.c *_embed.*

# End of file 
