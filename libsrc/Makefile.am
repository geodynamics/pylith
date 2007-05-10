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
	faults \
	materials \
	meshio \
	topology

lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	faults/BruneSlipFn.cc \
	faults/CohesiveTopology.cc \
	faults/EqKinSrc.cc \
	faults/Fault.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveKin.cc \
	faults/FaultCohesiveDyn.cc \
	faults/SlipTimeFn.cc \
	feassemble/ExplicitElasticity.cc \
	feassemble/ImplicitElasticity.cc \
	feassemble/Integrator.cc \
	feassemble/IntegratorElasticity.cc \
	feassemble/IntegratorExplicit.cc \
	feassemble/IntegratorImplicit.cc \
	feassemble/Quadrature.cc \
	feassemble/Quadrature1D.cc \
	feassemble/Quadrature1Din2D.cc \
	feassemble/Quadrature1Din3D.cc \
	feassemble/Quadrature2D.cc \
	feassemble/Quadrature2Din3D.cc \
	feassemble/Quadrature3D.cc \
	feassemble/ParameterManager.cc \
	materials/Material.cc \
	materials/ElasticStress1D.cc \
	materials/ElasticStrain1D.cc \
	materials/ElasticIsotropic3D.cc \
	materials/ElasticMaterial.cc \
	materials/ElasticPlaneStrain.cc \
	materials/ElasticPlaneStress.cc \
	meshio/BinaryIO.cc \
	meshio/GMVFile.cc \
	meshio/GMVFileAscii.cc \
	meshio/GMVFileBinary.cc \
	meshio/MeshIO.cc \
	meshio/MeshIOAscii.cc \
	meshio/MeshIOLagrit.cc \
	meshio/PsetFile.cc \
	meshio/PsetFileAscii.cc \
	meshio/PsetFileBinary.cc

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/MeshIOCubit.cc
endif

INCLUDES = -I$(top_builddir)/include
INCLUDES += $(PETSC_INCLUDE)

AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_SIEVE_FLAGS)

# End of file 
