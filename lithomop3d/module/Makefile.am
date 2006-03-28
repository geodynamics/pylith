## Process this file with automake to produce Makefile.in
# $Id: Makefile.am,v 1.1 2005/08/24 22:45:51 leif Exp $

noinst_LIBRARIES =
pkgpyexec_LTLIBRARIES =

if COND_EMBEDDING
    # static library
    noinst_LIBRARIES += liblithomop3dmodule.a
else
    # extension module (libtool)
    pkgpyexec_LTLIBRARIES += lithomop3dmodule.la
endif

# static library
liblithomop3dmodule_a_CXXFLAGS = $(AM_CXXFLAGS) # hack for automake
liblithomop3dmodule_a_SOURCES = $(sources)

# extension module (libtool)
lithomop3dmodule_la_SOURCES = $(sources)
lithomop3dmodule_la_LDFLAGS = -module
lithomop3dmodule_la_LIBADD = \
	$(top_builddir)/lithomop3d/liblithomop3d/liblithomop3d.la \
	-L$(PETSC_LIB_DIR) \
	$(PETSC_FORTRAN_LIB_BASIC) \
	$(PETSC_LIB_BASIC) \
	$(PETSC_EXTERNAL_LIB_BASIC) \
	$(FCLIBS) \
	-ljournal

CXX = $(MPICXX)
CXXLD = @CXX@
INCLUDES = $(PETSC_INCLUDE) -I$(PETSC_DIR)/src/dm/mesh/sieve -I$(PYTHON_INCDIR) $(MPIINCLUDES)
sources = \
	array.cc \
	array.h \
	autoprestr.cc \
	autoprestr.h \
	bindings.cc \
	bindings.h \
	elastc.cc \
	elastc.h \
	exceptionhandler.cc \
	exceptionhandler.h \
	exceptions.cc \
	exceptions.h \
	lithomop3dmodule.cc \
	lithomop3dmodule.h \
	lithomop3d_externs.h \
	misc.cc \
	misc.h \
	numbering.cc \
	numbering.h \
	parser.cc \
	parser.h \
	scanner.cc \
	scanner.h \
	setup.cc \
	setup.h \
	sorting.cc \
	sorting.h \
	sparse.cc \
	sparse.h \
	viscos.cc \
	viscos.h \
	write_modelinfo.cc \
	write_modelinfo.h

## end of Makefile.am
