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
	$(PYTHON_EGG_LDFLAGS) \
	$(PETSC_LIB) \
	$(FCLIBS)

INCLUDES = $(PETSC_INCLUDE) -I$(PYTHON_INCDIR)
AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) $(PETSC_SIEVE_FLAGS)
sources = \
	array.cc \
	array.h \
	bindings.cc \
	bindings.h \
	exceptionhandler.cc \
	exceptionhandler.h \
	exceptions.cc \
	exceptions.h \
	libpylith3d.cc \
	libpylith3d.h \
	pylith3dmodule.cc \
	pylith3dmodule.h \
	pylith3d_externs.h \
	mesh.cc \
	mesh.h \
	misc.cc \
	misc.h \
	setup.cc \
	setup.h

if ENABLE_GREENSFNS
sources += interpolate.cc \
	interpolate.h
CPPFLAGS += -DGREENSFNS
endif

## end of Makefile.am
