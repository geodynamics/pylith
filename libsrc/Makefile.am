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
	bc/DirichletBC.cc \
	bc/DirichletBoundary.cc \
	bc/Neumann.cc \
	feassemble/CellGeometry.cc \
	feassemble/Constraint.cc \
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
	feassemble/QuadratureRefCell.cc \
	feassemble/QuadratureEngine.cc \
	feassemble/Quadrature0D.cc \
	feassemble/Quadrature1D.cc \
	feassemble/Quadrature1Din2D.cc \
	feassemble/Quadrature1Din3D.cc \
	feassemble/Quadrature2D.cc \
	feassemble/Quadrature2Din3D.cc \
	feassemble/Quadrature3D.cc \
	meshio/BinaryIO.cc \
	meshio/GMVFile.cc \
	meshio/GMVFileAscii.cc \
	meshio/GMVFileBinary.cc \
	meshio/MeshBuilder.cc \
	meshio/MeshIO.cc \
	meshio/MeshIOAscii.cc \
	meshio/MeshIOLagrit.cc \
	meshio/PsetFile.cc \
	meshio/PsetFileAscii.cc \
	meshio/PsetFileBinary.cc \
	topology/Mesh.cc \
	topology/MeshOps.cc \
	topology/SubMesh.cc \
	topology/SolutionFields.cc \
	utils/EventLogger.cc

# 	faults/BruneSlipFn.cc \
# 	faults/ConstRateSlipFn.cc \
# 	faults/CohesiveTopology.cc \
# 	faults/EqKinSrc.cc \
# 	faults/Fault.cc \
# 	faults/FaultCohesive.cc \
# 	faults/FaultCohesiveKin.cc \
# 	faults/FaultCohesiveDyn.cc \
# 	faults/LiuCosSlipFn.cc \
# 	faults/SlipTimeFn.cc \
# 	faults/StepSlipFn.cc \
# 	feassemble/ElasticityExplicit.cc \
# 	feassemble/ElasticityImplicit.cc \
# 	feassemble/IntegratorElasticity.cc \
# 	feassemble/Quadrature.cc \
# 	feassemble/Quadrature0D.cc \
# 	feassemble/Quadrature1D.cc \
# 	feassemble/Quadrature1Din2D.cc \
# 	feassemble/Quadrature1Din3D.cc \
# 	feassemble/Quadrature2D.cc \
# 	feassemble/Quadrature2Din3D.cc \
# 	feassemble/Quadrature3D.cc \
# 	materials/Material.cc \
# 	materials/ElasticStress1D.cc \
# 	materials/ElasticStrain1D.cc \
# 	materials/ElasticIsotropic3D.cc \
# 	materials/ElasticMaterial.cc \
# 	materials/ElasticPlaneStrain.cc \
# 	materials/ElasticPlaneStress.cc \
# 	materials/MaxwellIsotropic3D.cc \
# 	materials/GenMaxwellIsotropic3D.cc \
# 	materials/ViscoelasticMaxwell.cc \
# 	meshio/CellFilter.cc \
# 	meshio/CellFilterAvg.cc \
# 	meshio/DataWriter.cc \
# 	meshio/DataWriterVTK.cc \
# 	meshio/OutputManager.cc \
# 	meshio/OutputSolnSubset.cc \
# 	meshio/VertexFilter.cc \
# 	meshio/VertexFilterVecNorm.cc \
# 	meshio/UCDFaultFile.cc \
# 	topology/Distributor.cc \
# 	topology/FieldsManager.cc \
# 	topology/MeshRefiner.cc 

#	topology/RefineUniform.cc \
#	bc/AbsorbingDampers.cc

libpylith_la_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
libpylith_la_LIBADD = \
	-lspatialdata \
	$(PETSC_LIB)
if NO_UNDEFINED
libpylith_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES = -I$(top_builddir)/include
INCLUDES += $(PETSC_INCLUDE)

AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_SIEVE_FLAGS)

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/MeshIOCubit.cc
  libpylith_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif

# End of file 
