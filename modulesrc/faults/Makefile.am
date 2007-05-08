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

subpackage = faults
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = faultsmodule.la

faultsmodule_la_LDFLAGS = -module

faultsmodule_la_SOURCES = faults.pyxe

nodist_faultsmodule_la_SOURCES = \
	faults.c faults_embed.cpp faults_embed.h

faultsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	$(PETSC_LIB)

if ENABLE_CUBIT
  faultsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

faults.pyx faults_embed.cpp  faults_embed.h: faults.pyxe
	pyrexembed faults.pyxe
faults.pyxe: $(srcdir)/faults.pyxe.src
	cp $(srcdir)/faults.pyxe.src $@
faults_embed.cpp: faults_embed.h
faults_embed.h: faults.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = faults.pyxe faults.pyx faults.c *_embed.*

# End of file 
