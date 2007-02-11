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

SUBDIRS = \
	utils \
	feassemble \
	meshio

# materials

lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	feassemble/ExplicitElasticity.cc \
	feassemble/Integrator.cc \
	feassemble/IntegratorExplicit.cc \
	feassemble/Quadrature.cc \
	feassemble/Quadrature1D.cc \
	feassemble/Quadrature1Din2D.cc \
	feassemble/Quadrature1Din3D.cc \
	feassemble/Quadrature2D.cc \
	feassemble/Quadrature2Din3D.cc \
	feassemble/Quadrature3D.cc \
	feassemble/ParameterManager.cc \
	meshio/MeshIO.cc \
	meshio/MeshIOAscii.cc


INCLUDES = -I$(top_builddir)/include
INCLUDES += $(PETSC_INCLUDE)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS)

# End of file 
