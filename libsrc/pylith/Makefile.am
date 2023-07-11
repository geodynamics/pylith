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
# Copyright (c) 2010-2022 University of California, Davis
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
	fekernels


lib_LTLIBRARIES = libpylith.la

libpylith_la_SOURCES = \
	bc/BoundaryCondition.cc \
	bc/TimeDependentAuxiliaryFactory.cc \
	bc/DirichletTimeDependent.cc \
	bc/DirichletUserFn.cc \
	bc/NeumannTimeDependent.cc \
	bc/NeumannUserFn.cc \
	bc/AbsorbingDampers.cc \
	bc/AbsorbingDampersAuxiliaryFactory.cc \
	faults/FaultCohesive.cc \
	faults/FaultCohesiveKin.cc \
	faults/FaultCohesiveImpulses.cc \
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
	faults/FaultOps.cc \
	feassemble/PhysicsImplementation.cc \
	feassemble/Integrator.cc \
	feassemble/IntegratorDomain.cc \
	feassemble/IntegratorBoundary.cc \
	feassemble/IntegratorInterface.cc \
	feassemble/IntegrationData.cc \
	feassemble/InterfacePatches.cc \
	feassemble/UpdateStateVars.cc \
	feassemble/JacobianValues.cc \
	feassemble/Constraint.cc \
	feassemble/ConstraintSpatialDB.cc \
	feassemble/ConstraintUserFn.cc \
	feassemble/ConstraintSimple.cc \
	feassemble/AuxiliaryFactory.cc \
	fekernels/Tensor.cc \
	fekernels/IsotropicLinearGenMaxwell.cc \
	fekernels/IsotropicPowerLaw.cc \
	materials/Material.cc \
	materials/Elasticity.cc \
	materials/RheologyElasticity.cc \
	materials/AuxiliaryFactoryElasticity.cc \
	materials/DerivedFactoryElasticity.cc \
	materials/DerivedFactoryPoroelasticity.cc \
	materials/IsotropicLinearElasticity.cc \
	materials/IsotropicLinearMaxwell.cc \
	materials/IsotropicLinearGenMaxwell.cc \
	materials/IsotropicPowerLaw.cc \
	materials/IncompressibleElasticity.cc \
	materials/RheologyIncompressibleElasticity.cc \
	materials/IsotropicLinearIncompElasticity.cc \
	materials/AuxiliaryFactoryElastic.cc \
	materials/AuxiliaryFactoryViscoelastic.cc \
	materials/Poroelasticity.cc \
	materials/RheologyPoroelasticity.cc \
	materials/IsotropicLinearPoroelasticity.cc \
	materials/AuxiliaryFactoryPoroelasticity.cc \
	materials/AuxiliaryFactoryPoroelastic.cc \
	materials/Query.cc \
	meshio/BinaryIO.cc \
	meshio/MeshBuilder.cc \
	meshio/MeshIO.cc \
	meshio/MeshIOAscii.cc \
	meshio/MeshIOPetsc.cc \
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
	problems/GreensFns.cc \
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
	problems/ProgressMonitorStep.cc \
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
	utils/PetscOptions.cc \
	utils/PylithVersion.cc \
	utils/PetscVersion.cc \
	utils/DependenciesVersion.cc \
	utils/TestArray.cc



libpylith_la_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LDFLAGS) $(PYTHON_LA_LDFLAGS)
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


# End of file
