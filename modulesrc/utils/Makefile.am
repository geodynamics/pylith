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

subpackage = utils
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = petscmodule.la

petscmodule_la_LDFLAGS = -module -no-undefined -avoid-version \
	$(PYTHON_LA_LDFLAGS)

petscmodule_la_SOURCES = petsc.pyxe

nodist_petscmodule_la_SOURCES = \
	petsc.c petsc_embed.cpp petsc_embed.h

petscmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  petscmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

petsc.pyx petsc_embed.cpp  petsc_embed.h: petsc.pyxe
	pyrexembed petsc.pyxe
petsc.pyxe: $(srcdir)/petsc.pyxe.src
	cp $(srcdir)/petsc.pyxe.src $@
petsc_embed.cpp: petsc_embed.h
petsc_embed.h: petsc.pyx

.pyx.c:
	pyrexc $< $(PYREX_INCLUDES)

CLEANFILES = petsc.pyxe petsc.pyx petsc.c *_embed.*

# End of file 
