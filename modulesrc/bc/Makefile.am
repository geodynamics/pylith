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

subpackage = bc
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = bcmodule.la

bcmodule_la_LDFLAGS = -module -no-undefined -avoid-version \
	$(PYTHON_LA_LDFLAGS)

bcmodule_la_SOURCES = bc.pyxe

nodist_bcmodule_la_SOURCES = \
	bc.c bc_embed.cpp bc_embed.h

bcmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  bcmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

bc.pyx bc_embed.cpp  bc_embed.h: bc.pyxe
	pyrexembed bc.pyxe
bc.pyxe: $(srcdir)/bc.pyxe.src
	cp $(srcdir)/bc.pyxe.src $@
bc_embed.cpp: bc_embed.h
bc_embed.h: bc.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = bc.pyxe bc.pyx bc.c *_embed.*

# End of file 
