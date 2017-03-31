# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

SUBDIRS = \
	utils \
	topology \
	meshio \
	materials \
	bc \
	feassemble \
	faults \
	friction \
	problems

lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	bc/BoundaryCondition.cc \
	bc/BoundaryConditionPoints.cc \
	bc/BCIntegratorSubMesh.cc \
	bc/TimeDependent.cc \
	bc/TimeDependentPoints.cc \
	bc/DirichletBC.cc \
	bc/DirichletBoundary.cc \
	bc/Neumann.cc \
	bc/AbsorbingDampers.cc \
	bc/PointForce.cc \
	faults/Fault.cc \
	faults/TopologyOps.cc \
	faults/CohesiveTopology.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveLagrange.cc \
	faults/FaultCohesiveKin.cc \
	faults/FaultCohesiveDyn.cc \
	faults/FaultCohesiveImpulses.cc \
	faults/FaultCohesiveTract.cc \
	faults/EqKinSrc.cc \
	faults/SlipTimeFn.cc \
	faults/StepSlipFn.cc \
	faults/ConstRateSlipFn.cc \
	faults/BruneSlipFn.cc \
	faults/TimeHistorySlipFn.cc \
	faults/LiuCosSlipFn.cc \
	faults/TractPerturbation.cc \
	feassemble/CellGeometry.cc \
	feassemble/Constraint.cc \
	feassemble/GeometryLine2D.cc \
	feassemble/GeometryLine3D.cc \
	feassemble/GeometryTri2D.cc \
	feassemble/GeometryTri3D.cc \
	feassemble/GeometryTet3D.cc \
	feassemble/GeometryQuad2D.cc \
	feassemble/GeometryQuad3D.cc \
	feassemble/GeometryHex3D.cc \
	feassemble/Quadrature.cc \
	feassemble/QuadratureRefCell.cc \
	feassemble/QuadratureEngine.cc \
	feassemble/Quadrature1Din2D.cc \
	feassemble/Quadrature1Din3D.cc \
	feassemble/Quadrature2D.cc \
	feassemble/Quadrature2Din3D.cc \
	feassemble/Quadrature3D.cc \
	feassemble/Integrator.cc \
	feassemble/IntegratorElasticity.cc \
	feassemble/ElasticityImplicit.cc \
	feassemble/ElasticityExplicit.cc \
	feassemble/ElasticityExplicitTri3.cc \
	feassemble/ElasticityExplicitTet4.cc \
	feassemble/IntegratorElasticityLgDeform.cc \
	feassemble/ElasticityImplicitLgDeform.cc \
	feassemble/ElasticityExplicitLgDeform.cc \
	friction/FrictionModel.cc \
	friction/StaticFriction.cc \
	friction/SlipWeakening.cc \
	friction/SlipWeakeningTime.cc \
	friction/SlipWeakeningTimeStable.cc \
	friction/RateStateAgeing.cc \
	friction/TimeWeakening.cc \
	materials/Metadata.cc \
	materials/Material.cc \
	materials/ElasticMaterial.cc \
	materials/ElasticPlaneStrain.cc \
	materials/ElasticPlaneStress.cc \
	materials/ElasticIsotropic3D.cc \
	materials/ViscoelasticMaxwell.cc \
	materials/GenMaxwellIsotropic3D.cc \
	materials/GenMaxwellPlaneStrain.cc \
	materials/GenMaxwellQpQsIsotropic3D.cc \
	materials/MaxwellIsotropic3D.cc \
	materials/MaxwellPlaneStrain.cc \
	materials/PowerLaw3D.cc \
	materials/PowerLawPlaneStrain.cc \
	materials/DruckerPrager3D.cc \
	materials/DruckerPragerPlaneStrain.cc \
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
	meshio/OutputSolnSubset.cc \
	meshio/OutputSolnPoints.cc \
	meshio/CellFilter.cc \
	meshio/CellFilterAvg.cc \
	meshio/VertexFilter.cc \
	meshio/VertexFilterVecNorm.cc \
	meshio/DataWriter.cc \
	meshio/DataWriterVTK.cc \
	meshio/OutputManager.cc \
	problems/Formulation.cc \
	problems/Explicit.cc \
	problems/Implicit.cc \
	problems/Solver.cc \
	problems/SolverLinear.cc \
	problems/SolverNonlinear.cc \
	problems/SolverLumped.cc \
	topology/FieldBase.cc \
	topology/Jacobian.cc \
	topology/Mesh.cc \
	topology/MeshOps.cc \
	topology/Field.cc \
	topology/Fields.cc \
	topology/SolutionFields.cc \
	topology/Distributor.cc \
	topology/ReverseCuthillMcKee.cc \
	topology/RefineUniform.cc \
	utils/EventLogger.cc \
	utils/PylithVersion.cc \
	utils/PetscVersion.cc \
	utils/DependenciesVersion.cc \
	utils/TestArray.cc


libpylith_la_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)
libpylith_la_LIBADD = \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS = -I$(top_srcdir)/libsrc
AM_CPPFLAGS += $(PETSC_CC_INCLUDES)

AM_CPPFLAGS += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_SIEVE_FLAGS)

if ENABLE_HDF5
  libpylith_la_SOURCES += \
	meshio/HDF5.cc \
	meshio/DataWriterHDF5.cc \
	meshio/DataWriterHDF5Ext.cc \
	meshio/Xdmf.cc
  libpylith_la_LIBADD += -lhdf5
endif

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/ExodusII.cc \
	meshio/MeshIOCubit.cc
  libpylith_la_LIBADD += -lnetcdf
endif


# End of file 
