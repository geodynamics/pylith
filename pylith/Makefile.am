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

nobase_pkgpyexec_PYTHON = \
	__init__.py \
	apps/__init__.py \
	apps/PyLithApp.py \
	apps/PetscApplication.py \
	apps/EqInfoApp.py \
	bc/__init__.py \
	bc/BoundaryCondition.py \
	bc/Dirichlet.py \
	bc/DirichletTimeDependent.py \
	bc/Neumann.py \
	bc/NeumannTimeDependent.py \
	bc/AuxFieldsTimeDependent.py \
	bc/AbsorbingDampers.py \
	bc/ZeroDB.py \
	faults/__init__.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveKin.py \
	faults/KinSrc.py \
	faults/KinSrcStep.py \
	faults/KinSrcConstRate.py \
	faults/SingleRupture.py \
	feassemble/__init__.py \
	feassemble/IntegratorPointwise.py \
	feassemble/ConstraintPointwise.py \
	feassemble/Observer.py \
	feassemble/SingleObserver.py \
	feassemble/ObservedComponent.py \
	friction/__init__.py \
	materials/__init__.py \
	materials/Homogeneous.py \
	materials/Material.py \
	materials/IsotropicLinearElasticityPlaneStrain.py \
	materials/IsotropicLinearMaxwellPlaneStrain.py \
	materials/IsotropicLinearIncompElasticityPlaneStrain.py \
	materials/AuxFieldsIsotropicLinearElasticity.py \
	materials/AuxFieldsIsotropicLinearMaxwell.py \
	meshio/__init__.py \
	meshio/DataWriter.py \
	meshio/DataWriterVTK.py \
	meshio/MeshIOObj.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOLagrit.py \
	meshio/PointsList.py \
	meshio/FieldFilter.py \
	meshio/FieldFilterNone.py \
	meshio/FieldFilterProject.py \
	meshio/OutputManager.py \
	meshio/OutputSoln.py \
	meshio/OutputMaterial.py \
	meshio/OutputTrigger.py \
	meshio/OutputTriggerStep.py \
	meshio/OutputTriggerTime.py \
	mpi/__init__.py \
	mpi/Communicator.py \
	problems/__init__.py \
	problems/Problem.py \
	problems/TimeDependent.py \
	problems/Solution.py \
	problems/SolutionSubfield.py \
	problems/SolnDisp.py \
	problems/SolnDispLagrange.py \
	problems/SolnDispVel.py \
	problems/SolnDispVelLagrange.py \
	problems/SolnDispPres.py \
	problems/SolnDispPresLagrange.py \
	problems/SubfieldDisplacement.py \
	problems/SubfieldVelocity.py \
	problems/SubfieldPressure.py \
	problems/SubfieldTemperature.py \
	problems/SubfieldLagrangeFault.py \
	topology/__init__.py \
	topology/Distributor.py \
	topology/Mesh.py \
	topology/Field.py \
	topology/Fields.py \
	topology/AuxSubfield.py \
	topology/MeshGenerator.py \
	topology/MeshImporter.py \
	topology/MeshImporterDist.py \
	topology/MeshRefiner.py \
	topology/RefineUniform.py \
	topology/ReverseCuthillMcKee.py \
	utils/__init__.py \
	utils/EmptyBin.py \
	utils/EventLogger.py \
	utils/NullComponent.py \
	utils/PetscComponent.py \
	utils/PetscManager.py \
	utils/CollectVersionInfo.py \
	utils/DumpParameters.py \
	utils/DumpParametersAscii.py \
	utils/DumpParametersJson.py \
	utils/profiling.py \
	utils/testarray.py \
	tests/__init__.py \
	tests/UnitTestApp.py \
	tests/FullTestApp.py \
	tests/PhysicalProperties.py \
	tests/Solution.py \
	tests/SolutionPoints.py \
	tests/StateVariables.py \
	tests/Fault.py

if ENABLE_CUBIT
  nobase_pkgpyexec_PYTHON += meshio/MeshIOCubit.py
endif

if ENABLE_HDF5
  nobase_pkgpyexec_PYTHON += \
	meshio/DataWriterHDF5.py \
	meshio/DataWriterHDF5Ext.py \
	meshio/Xdmf.py
endif


# End of file
