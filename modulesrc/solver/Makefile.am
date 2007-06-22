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

subpackage = solver
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = solvermodule.la

solvermodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

solvermodule_la_SOURCES = solver.pyxe

nodist_solvermodule_la_SOURCES = \
	solver.c solver_embed.cpp solver_embed.h

solvermodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  solvermodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

solver.pyx solver_embed.cpp  solver_embed.h: solver.pyxe
	pyrexembed solver.pyxe
solver.pyxe: $(srcdir)/solver.pyxe.src
	cp $(srcdir)/solver.pyxe.src $@
solver_embed.cpp: solver_embed.h
solver_embed.h: solver.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = solver.pyxe solver.pyx solver.c *_embed.*

# End of file 
