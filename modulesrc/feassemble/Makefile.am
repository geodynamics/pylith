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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = feassemblemodule.la

feassemblemodule_la_LDFLAGS = -module -no-undefined -avoid-version \
	$(PYTHON_LA_LDFLAGS)

feassemblemodule_la_SOURCES = feassemble.pyxe

nodist_feassemblemodule_la_SOURCES = \
	feassemble.c feassemble_embed.cpp feassemble_embed.h

feassemblemodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  feassemblemodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

feassemble.pyx feassemble_embed.cpp  feassemble_embed.h: feassemble.pyxe
	pyrexembed feassemble.pyxe
feassemble.pyxe: $(srcdir)/feassemble.pyxe.src
	cp $(srcdir)/feassemble.pyxe.src $@
feassemble_embed.cpp: feassemble_embed.h
feassemble_embed.h: feassemble.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = feassemble.pyxe feassemble.pyx feassemble.c *_embed.*

# End of file 
