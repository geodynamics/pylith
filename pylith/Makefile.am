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
# Copyright (c) 2010-2016 University of California, Davis
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
	bc/AbsorbingDampers.py \
	bc/BoundaryCondition.py \
	bc/DirichletBC.py \
	bc/DirichletBoundary.py \
	bc/Neumann.py \
	bc/PointForce.py \
	bc/TimeDependent.py \
	bc/TimeDependentPoints.py \
	bc/ZeroDispDB.py \
	faults/__init__.py \
	faults/BruneSlipFn.py \
	faults/ConstRateSlipFn.py \
	faults/EqKinSrc.py \
	faults/LiuCosSlipFn.py \
	faults/SingleRupture.py \
	faults/SlipTimeFn.py \
	faults/StepSlipFn.py \
	faults/TimeHistorySlipFn.py \
	faults/TractPerturbation.py \
	faults/Fault.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveKin.py \
	faults/FaultCohesiveDyn.py \
	faults/FaultCohesiveImpulses.py \
	faults/FaultCohesiveTract.py \
	feassemble/__init__.py \
	feassemble/Constraint.py \
	feassemble/ElasticityExplicit.py \
	feassemble/ElasticityExplicitTet4.py \
	feassemble/ElasticityExplicitTri3.py \
	feassemble/ElasticityExplicitLgDeform.py \
	feassemble/ElasticityImplicit.py \
	feassemble/ElasticityImplicitLgDeform.py \
	feassemble/FIATQuadrature.py \
	feassemble/FIATLagrange.py \
	feassemble/FIATSimplex.py \
	feassemble/Integrator.py \
	feassemble/IntegratorElasticity.py \
	feassemble/IntegratorElasticityLgDeform.py \
	feassemble/ReferenceCell.py \
	feassemble/Quadrature.py \
	feassemble/CellGeometry.py \
	feassemble/IntegratorPointwise.py \
	feassemble/ConstraintPointwise.py \
	friction/__init__.py \
	friction/FrictionModel.py \
	friction/StaticFriction.py \
	friction/SlipWeakening.py \
	friction/SlipWeakeningTime.py \
	friction/SlipWeakeningTimeStable.py \
	friction/RateStateAgeing.py \
	friction/TimeWeakening.py \
	materials/__init__.py \
	materials/ElasticMaterial.py \
	materials/ElasticIsotropic3D.py \
	materials/ElasticPlaneStrain.py \
	materials/ElasticPlaneStress.py \
	materials/GenMaxwellIsotropic3D.py \
	materials/GenMaxwellPlaneStrain.py \
	materials/GenMaxwellQpQsIsotropic3D.py \
	materials/Homogeneous.py \
	materials/Material.py \
	materials/MaxwellIsotropic3D.py \
	materials/MaxwellPlaneStrain.py \
	materials/PowerLaw3D.py \
	materials/PowerLawPlaneStrain.py \
	materials/DruckerPrager3D.py \
	materials/DruckerPragerPlaneStrain.py \
	materials/MaterialNew.py \
	materials/IsotropicLinearElasticityPlaneStrain.py \
	materials/AuxFieldsIsotropicLinearElasticity.py \
	meshio/__init__.py \
	meshio/CellFilter.py \
	meshio/CellFilterAvg.py \
	meshio/DataWriter.py \
	meshio/DataWriterVTK.py \
	meshio/MeshIOObj.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOLagrit.py \
	meshio/OutputDirichlet.py \
	meshio/OutputManager.py \
	meshio/OutputSoln.py \
	meshio/OutputFaultKin.py \
	meshio/OutputFaultDyn.py \
	meshio/OutputFaultImpulses.py \
	meshio/OutputMatElastic.py \
	meshio/OutputNeumann.py \
	meshio/OutputSolnSubset.py \
	meshio/OutputSolnPoints.py \
	meshio/PointsList.py \
	meshio/SingleOutput.py \
	meshio/VertexFilter.py \
	meshio/VertexFilterVecNorm.py \
	mpi/__init__.py \
	mpi/Communicator.py \
	perf/__init__.py \
	perf/Memory.py \
	perf/Logger.py \
	perf/MemoryLogger.py \
	perf/Mesh.py \
	perf/Fault.py \
	perf/Material.py \
	perf/VertexGroup.py \
	perf/Field.py \
	perf/GlobalOrder.py \
	perf/Jacobian.py \
	problems/__init__.py \
	problems/Explicit.py \
	problems/ExplicitTri3.py \
	problems/ExplicitTet4.py \
	problems/ExplicitLgDeform.py \
	problems/Formulation.py \
	problems/Implicit.py \
	problems/ImplicitLgDeform.py \
	problems/Problem.py \
	problems/ProblemNew.py \
	problems/Solver.py \
	problems/SolverLinear.py \
	problems/SolverNonlinear.py \
	problems/SolverLumped.py \
	problems/TimeDependent.py \
	problems/TimeDependentNew.py \
	problems/GreensFns.py \
	problems/TimeStep.py \
	problems/TimeStepAdapt.py \
	problems/TimeStepUniform.py \
	problems/TimeStepUser.py \
	problems/ProgressMonitor.py \
	problems/ProgressMonitorStep.py \
	problems/ProgressMonitorTime.py \
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
	topology/SolutionFields.py \
	topology/AuxSubfield.py \
	topology/Jacobian.py \
	topology/JacobianViewer.py \
	topology/MeshGenerator.py \
	topology/MeshImporter.py \
	topology/MeshImporterDist.py \
	topology/MeshRefiner.py \
	topology/RefineUniform.py \
	topology/ReverseCuthillMcKee.py \
	utils/__init__.py \
	utils/CheckpointTimer.py \
	utils/CppData.py \
	utils/EmptyBin.py \
	utils/EventLogger.py \
	utils/NullComponent.py \
	utils/PetscComponent.py \
	utils/PetscManager.py \
	utils/VTKDataReader.py \
	utils/importing.py \
	utils/profiling.py \
	utils/testarray.py \
	tests/__init__.py \
	tests/CppTestData.py \
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
