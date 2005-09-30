## Process this file with automake to produce Makefile.in
# $Id: Makefile.am,v 1.1 2005/08/24 22:45:51 leif Exp $

CXX = $(MPICXX)
CXXLD = @CXX@
INCLUDES = $(PETSC_INCLUDE) -I$(PYTHON_INCDIR) $(MPIINCLUDES)
pkgpyexec_LTLIBRARIES = lithomop3dmodule.la
lithomop3dmodule_la_SOURCES = \
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
	lithomop3d_externs.h \
	lithomop3dmodule.cc \
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
lithomop3dmodule_la_LDFLAGS = -module
lithomop3dmodule_la_LIBADD = \
	$(top_builddir)/lithomop3d/liblithomop3d/liblithomop3d.la \
	-L$(PETSC_LIB_DIR) \
	$(PETSC_FORTRAN_LIB_BASIC) \
	$(PETSC_LIB_BASIC) \
	$(PETSC_EXTERNAL_LIB_BASIC) \
	-ljournal

## end of Makefile.am
