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
	topology \
	meshio \
	feassemble \
	bc \
	faults \
	materials

lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	bc/BoundaryCondition.cc \
	bc/Dirichlet.cc \
	faults/BruneSlipFn.cc \
	faults/CohesiveTopology.cc \
	faults/EqKinSrc.cc \
	faults/Fault.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveKin.cc \
	faults/FaultCohesiveDyn.cc \
	faults/SlipTimeFn.cc \
	feassemble/CellGeometry.cc \
	feassemble/Constraint.cc \
	feassemble/Elasticity.cc \
	feassemble/ElasticityExplicit.cc \
	feassemble/ElasticityImplicit.cc \
	feassemble/GeometryPoint1D.cc \
	feassemble/GeometryPoint2D.cc \
	feassemble/GeometryPoint3D.cc \
	feassemble/GeometryLine1D.cc \
	feassemble/GeometryLine2D.cc \
	feassemble/GeometryLine3D.cc \
	feassemble/GeometryTri2D.cc \
	feassemble/GeometryTri3D.cc \
	feassemble/GeometryTet3D.cc \
	feassemble/GeometryQuad2D.cc \
	feassemble/GeometryQuad3D.cc \
	feassemble/GeometryHex3D.cc \
	feassemble/Integrator.cc \
	feassemble/Quadrature.cc \
	feassemble/Quadrature0D.cc \
	feassemble/Quadrature1D.cc \
	feassemble/Quadrature1Din2D.cc \
	feassemble/Quadrature1Din3D.cc \
	feassemble/Quadrature2D.cc \
	feassemble/Quadrature2Din3D.cc \
	feassemble/Quadrature3D.cc \
	materials/Material.cc \
	materials/ElasticStress1D.cc \
	materials/ElasticStrain1D.cc \
	materials/ElasticIsotropic3D.cc \
	materials/ElasticMaterial.cc \
	materials/ElasticPlaneStrain.cc \
	materials/ElasticPlaneStress.cc \
	materials/MaxwellIsotropic3D.cc \
	meshio/BinaryIO.cc \
	meshio/GMVFile.cc \
	meshio/GMVFileAscii.cc \
	meshio/GMVFileBinary.cc \
	meshio/MeshIO.cc \
	meshio/MeshIOAscii.cc \
	meshio/MeshIOLagrit.cc \
	meshio/PsetFile.cc \
	meshio/PsetFileAscii.cc \
	meshio/PsetFileBinary.cc \
	meshio/SolutionIO.cc \
	meshio/SolutionIOVTK.cc \
	topology/FieldsManager.cc

libpylith_la_LDFLAGS = -no-undefined $(PYTHON_LA_LDFLAGS)
libpylith_la_LIBADD = \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

INCLUDES = -I$(top_builddir)/include
INCLUDES += $(PETSC_INCLUDE)

AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_SIEVE_FLAGS)

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/MeshIOCubit.cc
  libpylith_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

# End of file 
