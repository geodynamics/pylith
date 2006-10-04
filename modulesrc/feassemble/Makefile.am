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

feassemblemodule_la_LDFLAGS = -module

feassemblemodule_la_SOURCES = feassemble.pyxe

nodist_feassemblemodule_la_SOURCES = \
	feassemble.c feassemble_embed.cpp feassemble_embed.h

feassemblemodule_la_LIBADD = \
	$(top_builddir)/libsrc/feassemble/libpylithfeassemble.la \
	$(PETSC_LIB)

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

feassemble.pyx feassemble_embed.cpp  feassemble_embed.h: feassemble.pyxe
	if [ "${VPATH}" != "" ]; then cp $< .; fi && pyrexembed feassemble.pyxe && if [ "${VPATH}" != "" ]; then rm -f feassemble.pyxe; fi
feassemble_embed.cpp: feassemble_embed.h
feassemble_embed.h: feassemble.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = feassemble.pyx feassemble.c *_embed.*

# End of file 
