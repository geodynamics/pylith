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

if ENABLE_CUBIT
make-manifest:
else
make-manifest:
	echo "exclude meshio/MeshIOCubit.py" > $@
endif

install-exec-local: make-manifest
	$(PYTHON) -m pip install $(top_srcdir)

EXTRA_DIST = \
	__init__.py \
	apps/ConfigSearchApp.py \
	apps/EqInfoApp.py \
	apps/PetscApplication.py \
	apps/PyLithApp.py \
	apps/RunnerApp.py \
	apps/__init__.py \
	bc/AbsorbingDampers.py \
	bc/AuxSubfieldsAbsorbingDampers.py \
	bc/AuxSubfieldsTimeDependent.py \
	bc/BoundaryCondition.py \
	bc/DirichletTimeDependent.py \
	bc/NeumannTimeDependent.py \
	bc/ZeroDB.py \
	bc/__init__.py \
	faults/AuxSubfieldsFault.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveImpulses.py \
	faults/FaultCohesiveKin.py \
	faults/KinSrc.py \
	faults/KinSrcBrune.py \
	faults/KinSrcConstRate.py \
	faults/KinSrcLiuCos.py \
	faults/KinSrcRamp.py \
	faults/KinSrcStep.py \
	faults/KinSrcTimeHistory.py \
	faults/SingleRupture.py \
	faults/__init__.py \
	feassemble/__init__.py \
	friction/__init__.py \
	materials/AuxSubfieldsElasticity.py \
	materials/AuxSubfieldsIsotropicLinearElasticity.py \
	materials/AuxSubfieldsIsotropicLinearGenMaxwell.py \
	materials/AuxSubfieldsIsotropicLinearMaxwell.py \
	materials/AuxSubfieldsIsotropicLinearPoroelasticity.py \
	materials/AuxSubfieldsIsotropicPowerLaw.py \
	materials/AuxSubfieldsPoroelasticity.py \
	materials/DerivedSubfieldsElasticity.py \
	materials/Elasticity.py \
	materials/Homogeneous.py \
	materials/IncompressibleElasticity.py \
	materials/IsotropicLinearElasticity.py \
	materials/IsotropicLinearGenMaxwell.py \
	materials/IsotropicLinearIncompElasticity.py \
	materials/IsotropicLinearMaxwell.py \
	materials/IsotropicLinearPoroelasticity.py \
	materials/IsotropicPowerLaw.py \
	materials/Material.py \
	materials/Poroelasticity.py \
	materials/RheologyElasticity.py \
	materials/RheologyIncompressibleElasticity.py \
	materials/RheologyPoroelasticity.py \
	materials/__init__.py \
	meshio/DataWriter.py \
	meshio/DataWriterHDF5.py \
	meshio/DataWriterHDF5Ext.py \
	meshio/DataWriterVTK.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOCubit.py \
	meshio/MeshIOObj.py \
	meshio/MeshIOPetsc.py \
	meshio/OutputObserver.py \
	meshio/OutputPhysics.py \
	meshio/OutputSoln.py \
	meshio/OutputSolnBoundary.py \
	meshio/OutputSolnDomain.py \
	meshio/OutputSolnPoints.py \
	meshio/OutputTrigger.py \
	meshio/OutputTriggerStep.py \
	meshio/OutputTriggerTime.py \
	meshio/PointsList.py \
	meshio/Xdmf.py \
	meshio/__init__.py \
	meshio/gmsh_utils.py \
	mpi/Communicator.py \
	mpi/__init__.py \
	problems/GreensFns.py \
	problems/InitialCondition.py \
	problems/InitialConditionDomain.py \
	problems/InitialConditionPatch.py \
	problems/Physics.py \
	problems/Problem.py \
	problems/ProblemDefaults.py \
	problems/ProgressMonitor.py \
	problems/ProgressMonitorStep.py \
	problems/ProgressMonitorTime.py \
	problems/SingleObserver.py \
	problems/SolnDisp.py \
	problems/SolnDispLagrange.py \
	problems/SolnDispPres.py \
	problems/SolnDispPresLagrange.py \
	problems/SolnDispPresTracStrain.py \
	problems/SolnDispPresTracStrainVelPdotTdot.py \
	problems/SolnDispPresVel.py \
	problems/SolnDispVel.py \
	problems/SolnDispVelLagrange.py \
	problems/Solution.py \
	problems/SolutionSubfield.py \
	problems/SubfieldDisplacement.py \
	problems/SubfieldLagrangeFault.py \
	problems/SubfieldPressure.py \
	problems/SubfieldPressureDot.py \
	problems/SubfieldTemperature.py \
	problems/SubfieldTraceStrain.py \
	problems/SubfieldTraceStrainDot.py \
	problems/SubfieldVelocity.py \
	problems/TimeDependent.py \
	problems/__init__.py \
	testing/FullTestApp.py \
	testing/SolutionPoints.py \
	testing/UnitTestApp.py \
	testing/__init__.py \
	topology/Distributor.py \
	topology/Field.py \
	topology/Mesh.py \
	topology/MeshGenerator.py \
	topology/MeshImporter.py \
	topology/MeshImporterDist.py \
	topology/MeshRefiner.py \
	topology/RefineUniform.py \
	topology/ReverseCuthillMcKee.py \
	topology/Subfield.py \
	topology/__init__.py \
	utils/CollectVersionInfo.py \
	utils/DumpParameters.py \
	utils/DumpParametersAscii.py \
	utils/DumpParametersJson.py \
	utils/EmptyBin.py \
	utils/EventLogger.py \
	utils/NullComponent.py \
	utils/PetscComponent.py \
	utils/PetscDefaults.py \
	utils/PetscManager.py \
	utils/PropertyList.py \
	utils/SimulationMetadata.py \
	utils/__init__.py \
	utils/converters.py \
	utils/profiling.py \
	utils/validators.py 


# End of file
