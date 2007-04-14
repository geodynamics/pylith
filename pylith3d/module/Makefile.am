## Process this file with automake to produce Makefile.in
# $Id: Makefile.am,v 1.1 2005/08/24 22:45:51 leif Exp $

noinst_LIBRARIES =
pkgpyexec_LTLIBRARIES =

if COND_EMBEDDING
    # static library
    noinst_LIBRARIES += libPyLithLibmodule.a libPyLithMeshLibmodule.a
else
    # extension module (libtool)
    pkgpyexec_LTLIBRARIES += PyLithLibmodule.la PyLithMeshLibmodule.la
endif

# static libraries

libPyLithLibmodule_a_CFLAGS = $(AM_CFLAGS) # hack for automake
libPyLithLibmodule_a_CXXFLAGS = $(AM_CXXFLAGS) # hack for automake
libPyLithLibmodule_a_SOURCES = $(sources)

libPyLithMeshLibmodule_a_CFLAGS = $(AM_CFLAGS) # hack for automake
libPyLithMeshLibmodule_a_CXXFLAGS = $(AM_CXXFLAGS) # hack for automake
libPyLithMeshLibmodule_a_SOURCES = $(meshsources)

# extension modules (libtool)

PyLithLibmodule_la_SOURCES = $(sources)
PyLithLibmodule_la_LDFLAGS = -module
PyLithLibmodule_la_LIBADD = \
	$(top_builddir)/pylith3d/libpylith3d/libpylith3d.la \
	$(PYTHON_EGG_LDFLAGS) \
	$(PETSC_LIB) \
	$(FCLIBS)

PyLithMeshLibmodule_la_SOURCES = $(meshsources)
PyLithMeshLibmodule_la_LDFLAGS = -module
PyLithMeshLibmodule_la_LIBADD = \
	$(PYTHON_EGG_LDFLAGS) \
	$(PETSC_LIB) \
	$(FCLIBS)

INCLUDES = $(PETSC_INCLUDE) -I$(PYTHON_INCDIR)
AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) $(PETSC_SIEVE_FLAGS)

sources = \
	array.pyx \
	constants.pxd \
	elementtypes.cc \
	elementtypes.h \
	elementtypes.pxd \
	exceptionhandler.cc \
	exceptionhandler.pxd \
	PyLithLib.cc \
	setup.pyx

meshsources = \
	mesh.cc \
	mesh.h \
	PyLithMeshLib.cc \
	PyLithMeshLib.pxd

if ENABLE_GREENSFNS
meshsources += interpolate.cc \
	interpolate.h
CPPFLAGS += -DGREENSFNS
endif


mesh.cc: PyLithMeshLib.h


PYREX_FLAGS = -I$(srcdir) -I$(top_srcdir)/pylith3d/libpylith3d
pyxincludes = \
	$(srcdir)/constants.pxd \
	$(srcdir)/elementtypes.pxd \
	$(top_srcdir)/pylith3d/libpylith3d/libpylith3d.pxd \
	$(srcdir)/petsc.pxd \
	$(srcdir)/PyLithMeshLib.pxd \
	$(srcdir)/array.pyx \
	$(srcdir)/exceptionhandler.pxd \
	$(srcdir)/setup.pyx


PyLithLib.cc: $(srcdir)/PyLithLib.pyx $(pyxincludes)
	pyrexc $(PYREX_FLAGS) $< -o PyLithLib.cc

PyLithMeshLib.cc: $(srcdir)/PyLithMeshLib.pyx $(srcdir)/petsc.pxd
	pyrexc $(PYREX_FLAGS) $< -o PyLithMeshLib.cc

PyLithMeshLib.h: PyLithMeshLib.pxd

EXTRA_DIST = \
	$(srcdir)/PyLithLib.pyx \
	$(srcdir)/PyLithMeshLib.pyx \
	PyLithLib.h \
	PyLithMeshLib.h


## end of Makefile.am
