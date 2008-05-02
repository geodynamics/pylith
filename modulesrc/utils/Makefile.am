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

subpkgpyexec_LTLIBRARIES = petscmodule.la utilsmodule.la

petscmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
dist_petscmodule_la_SOURCES = petsc.pyxe.src
nodist_petscmodule_la_SOURCES = \
	petsc.pyxe \
	petsc.c petsc_embed.cpp petsc_embed.h
petscmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB)
if ENABLE_CUBIT
  petscmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
  petscmodule_la_LIBADD += $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif


utilsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
dist_utilsmodule_la_SOURCES = utils.pyxe.src
nodist_utilsmodule_la_SOURCES = \
	utils.pyxe \
	utils.c utils_embed.cpp utils_embed.h
utilsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB)
if ENABLE_CUBIT
  utilsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
  utilsmodule_la_LIBADD += $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif



INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

petsc.pyx petsc_embed.cpp  petsc_embed.h: petsc.pyxe
	pyrexembed petsc.pyxe
petsc.pyxe: $(srcdir)/petsc.pyxe.src
	cp $(srcdir)/petsc.pyxe.src $@
petsc_embed.cpp: petsc_embed.h
petsc_embed.h: petsc.pyx

utils.pyx utils_embed.cpp  utils_embed.h: utils.pyxe
	pyrexembed utils.pyxe
utils.pyxe: $(srcdir)/utils.pyxe.src
	cp $(srcdir)/utils.pyxe.src $@
utils_embed.cpp: utils_embed.h
utils_embed.h: utils.pyx

.pyx.c:
	pyrexc $< $(PYREX_INCLUDES)

CLEANFILES = *_embed.* \
	petsc.pyxe petsc.pyx petsc.c \
	utils.pyxe utils.pyx utils.c

# End of file 
