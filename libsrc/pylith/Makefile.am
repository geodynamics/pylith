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
	faults \
	friction \
	feassemble \
	fekernels \
	problems


lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	bc/BoundaryCondition.cc \
	bc/TimeDependentAuxiliaryFactory.cc \
	bc/ConstraintBoundary.cc \
	bc/DirichletTimeDependent.cc \
	bc/IntegratorBoundary.cc \
	bc/NeumannTimeDependent.cc \
	bc/AbsorbingDampers.cc \
	bc/AbsorbingDampersAuxiliaryFactory.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveStub.cc \
	faults/FaultCohesiveKin.cc \
	faults/AuxiliaryFactory.cc \
	faults/KinSrc.cc \
	faults/KinSrcStep.cc \
	faults/KinSrcConstRate.cc \
	faults/KinSrcAuxiliaryFactory.cc \
	faults/TopologyOps.cc \
	feassemble/IntegratorPointwise.cc \
	feassemble/ConstraintPointwise.cc \
	feassemble/Observer.cc \
	feassemble/ObservedComponent.cc \
	feassemble/AuxiliaryFactory.cc \
	fekernels/Elasticity.cc \
	fekernels/ElasticityPlaneStrain.cc \
	fekernels/IsotropicLinearElasticityPlaneStrain.cc \
	fekernels/Viscoelastic.cc \
	fekernels/IsotropicLinearMaxwellPlaneStrain.cc \
	fekernels/IsotropicLinearGenMaxwellPlaneStrain.cc \
	fekernels/IncompressibleElasticity.cc \
	fekernels/IsotropicLinearIncompElasticityPlaneStrain.cc \
	fekernels/DispVel.cc \
	fekernels/TimeDependentFn.cc \
	fekernels/NeumannTimeDependent.cc \
	fekernels/AbsorbingDampers.cc \
	fekernels/FaultCohesiveKin.cc \
	materials/Material.cc \
	materials/AuxiliaryFactory.cc \
	materials/IsotropicLinearElasticityPlaneStrain.cc \
	materials/IsotropicLinearMaxwellPlaneStrain.cc \
	materials/IsotropicLinearGenMaxwellPlaneStrain.cc \
	materials/IsotropicLinearIncompElasticityPlaneStrain.cc \
	materials/Query.cc \
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
	meshio/OutputSolnBoundary.cc \
	meshio/OutputSolnPoints.cc \
	meshio/FieldFilter.cc \
	meshio/FieldFilterNone.cc \
	meshio/FieldFilterProject.cc \
	meshio/DataWriter.cc \
	meshio/DataWriterVTK.cc \
	meshio/OutputManager.cc \
	meshio/OutputIntegrator.cc \
	meshio/OutputConstraint.cc \
	meshio/OutputSoln.cc \
	meshio/OutputTrigger.cc \
	meshio/OutputTriggerStep.cc \
	meshio/OutputTriggerTime.cc \
	problems/Problem.cc \
	problems/TimeDependent.cc \
	problems/SolutionFactory.cc \
	topology/FieldBase.cc \
	topology/Mesh.cc \
	topology/MeshOps.cc \
	topology/Field.cc \
	topology/FieldOps.cc \
	topology/FieldQuery.cc \
	topology/Fields.cc \
	topology/Distributor.cc \
	topology/ReverseCuthillMcKee.cc \
	topology/RefineUniform.cc \
	utils/EventLogger.cc \
	utils/PyreComponent.cc \
	utils/GenericComponent.cc \
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
	meshio/Xdmf.cc \
	meshio/DataWriterHDF5.cc \
	meshio/DataWriterHDF5Ext.cc
  libpylith_la_LIBADD += -lhdf5
endif

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/ExodusII.cc \
	meshio/MeshIOCubit.cc
  libpylith_la_LIBADD += -lnetcdf
endif


# End of file
