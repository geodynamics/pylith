## Process this file with automake to produce Makefile.in
# $Id: Makefile.am,v 1.1 2005/08/24 22:45:51 leif Exp $

noinst_LIBRARIES =
pkgpyexec_LTLIBRARIES =

if COND_EMBEDDING
    # static library
    noinst_LIBRARIES += libpylith3dmodule.a
else
    # extension module (libtool)
    pkgpyexec_LTLIBRARIES += pylith3dmodule.la
endif

# static library
libpylith3dmodule_a_CXXFLAGS = $(AM_CXXFLAGS) # hack for automake
libpylith3dmodule_a_SOURCES = $(sources)

# extension module (libtool)
pylith3dmodule_la_SOURCES = $(sources)
pylith3dmodule_la_LDFLAGS = -module
pylith3dmodule_la_LIBADD = \
	$(top_builddir)/pylith3d/libpylith3d/libpylith3d.la \
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
	pylith3dmodule.cc \
	pylith3dmodule.h \
	pylith3d_externs.h \
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
