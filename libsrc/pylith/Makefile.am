# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

SUBDIRS = \
	utils \
	topology \
	feassemble \
	problems \
	meshio \
	materials \
	bc \
	faults \
	friction \
	fekernels \
	sources \
	testing


lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	bc/BoundaryCondition.cc \
	bc/TimeDependentAuxiliaryFactory.cc \
	bc/DirichletTimeDependent.cc \
	bc/DirichletUserFn.cc \
	bc/NeumannTimeDependent.cc \
	bc/AbsorbingDampers.cc \
	bc/AbsorbingDampersAuxiliaryFactory.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveKin.cc \
	faults/AuxiliaryFactoryKinematic.cc \
	faults/KinSrc.cc \
	faults/KinSrcStep.cc \
	faults/KinSrcRamp.cc \
	faults/KinSrcConstRate.cc \
	faults/KinSrcBrune.cc \
	faults/KinSrcLiuCos.cc \
	faults/KinSrcTimeHistory.cc \
	faults/KinSrcAuxiliaryFactory.cc \
	faults/TopologyOps.cc \
	feassemble/PhysicsImplementation.cc \
	feassemble/Integrator.cc \
	feassemble/IntegratorDomain.cc \
	feassemble/IntegratorBoundary.cc \
	feassemble/IntegratorInterface.cc \
	feassemble/InterfacePatches.cc \
	feassemble/UpdateStateVars.cc \
	feassemble/Constraint.cc \
	feassemble/ConstraintSpatialDB.cc \
	feassemble/ConstraintUserFn.cc \
	feassemble/ConstraintSimple.cc \
	feassemble/AuxiliaryFactory.cc \
	fekernels/Elasticity.cc \
	fekernels/IsotropicLinearElasticity.cc \
	fekernels/Viscoelasticity.cc \
	fekernels/IsotropicLinearMaxwell.cc \
	fekernels/IsotropicLinearGenMaxwell.cc \
	fekernels/IsotropicPowerLaw.cc \
	fekernels/IncompressibleElasticity.cc \
	fekernels/IsotropicLinearIncompElasticity.cc \
	fekernels/Poroelasticity.cc \
	fekernels/IsotropicLinearPoroelasticity.cc \
	fekernels/Solution.cc \
	fekernels/DispVel.cc \
	fekernels/TimeDependentFn.cc \
	fekernels/NeumannTimeDependent.cc \
	fekernels/AbsorbingDampers.cc \
	fekernels/FaultCohesiveKin.cc \
	fekernels/WellboreSource.cc \
	materials/Material.cc \
	materials/Elasticity.cc \
	materials/RheologyElasticity.cc \
	materials/AuxiliaryFactoryElasticity.cc \
	materials/DerivedFactoryElasticity.cc \
	materials/IsotropicLinearElasticity.cc \
	materials/IsotropicLinearMaxwell.cc \
	materials/IsotropicLinearGenMaxwell.cc \
	materials/IsotropicPowerLaw.cc \
	materials/AuxiliaryFactoryElastic.cc \
	materials/AuxiliaryFactoryViscoelastic.cc \
	materials/AuxiliaryFactoryPoroelasticity.cc \
	materials/AuxiliaryFactoryPoroelastic.cc \
	materials/IncompressibleElasticity.cc \
	materials/RheologyIncompressibleElasticity.cc \
	materials/IsotropicLinearIncompElasticity.cc \
	materials/Poroelasticity.cc \
	materials/RheologyPoroelasticity.cc \
	materials/IsotropicLinearPoroelasticity.cc \
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
	meshio/DataWriter.cc \
	meshio/HDF5.cc \
	meshio/Xdmf.cc \
	meshio/DataWriterHDF5.cc \
	meshio/DataWriterHDF5Ext.cc \
	meshio/DataWriterVTK.cc \
	meshio/OutputObserver.cc \
	meshio/OutputSubfield.cc \
	meshio/OutputSoln.cc \
	meshio/OutputSolnDomain.cc \
	meshio/OutputSolnBoundary.cc \
	meshio/OutputSolnPoints.cc \
	meshio/OutputPhysics.cc \
	meshio/OutputTrigger.cc \
	meshio/OutputTriggerStep.cc \
	meshio/OutputTriggerTime.cc \
	problems/Problem.cc \
	problems/TimeDependent.cc \
	problems/SolutionFactory.cc \
	problems/ObserverSoln.cc \
	problems/ObserversSoln.cc \
	problems/Physics.cc \
	problems/ObserverPhysics.cc \
	problems/ObserversPhysics.cc \
	problems/InitialCondition.cc \
	problems/InitialConditionDomain.cc \
	problems/InitialConditionPatch.cc \
	problems/ProgressMonitor.cc \
	problems/ProgressMonitorTime.cc \
	sources/Source.cc \
	sources/WellboreSource.cc \
	sources/AuxiliaryFactoryWellboreSource.cc \
	topology/Mesh.cc \
	topology/MeshOps.cc \
	topology/FieldBase.cc \
	topology/Field.cc \
	topology/FieldFactory.cc \
	topology/FieldOps.cc \
	topology/FieldQuery.cc \
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
	-lhdf5 \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS = -I$(top_srcdir)/libsrc
AM_CPPFLAGS += $(PETSC_CC_INCLUDES)

AM_CPPFLAGS += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR)

if ENABLE_CUBIT
  libpylith_la_SOURCES += \
	meshio/ExodusII.cc \
	meshio/MeshIOCubit.cc
  libpylith_la_LIBADD += -lnetcdf
endif

if ENABLE_TESTING
libpylith_la_SOURCES += \
	testing/StubMethodTracker.cc \
	testing/FieldTester.cc \
	testing/PhysicsStub.cc \
	testing/PhysicsImplementationStub.cc \
	testing/ProgressMonitorStub.cc \
	testing/ObserverPhysicsStub.cc \
	testing/ObserverSolnStub.cc \
	testing/FaultCohesiveStub.cc \
	testing/MMSTest.cc \
	testing/TestDriver.cc

libpylith_la_LIBADD += -lcppunit
endif



# End of file
