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

nobase_pkgpyexec_PYTHON = \
	PyLithApp.py \
	__init__.py \
	bc/__init__.py \
	bc/AbsorbingDampers.py \
	bc/BoundaryCondition.py \
	bc/DirichletBC.py \
	bc/DirichletBoundary.py \
	bc/FixedDOFDB.py \
	bc/Neumann.py \
	faults/__init__.py \
	faults/BruneSlipFn.py \
	faults/ConstRateSlipFn.py \
	faults/EqKinSrc.py \
	faults/Fault.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveKin.py \
	faults/LiuCosSlipFn.py \
	faults/SingleRupture.py \
	faults/SlipTimeFn.py \
	faults/StepSlipFn.py \
	feassemble/__init__.py \
	feassemble/Constraint.py \
	feassemble/ElasticityExplicit.py \
	feassemble/ElasticityImplicit.py \
	feassemble/FIATLagrange.py \
	feassemble/FIATSimplex.py \
	feassemble/Integrator.py \
	feassemble/IntegratorElasticity.py \
	feassemble/ReferenceCell.py \
	feassemble/Quadrature.py \
	feassemble/CellGeometry.py \
	materials/__init__.py \
	materials/ElasticMaterial.py \
	materials/ElasticIsotropic3D.py \
	materials/ElasticPlaneStrain.py \
	materials/ElasticPlaneStress.py \
	materials/ElasticStrain1D.py \
	materials/ElasticStress1D.py \
	materials/GenMaxwellIsotropic3D.py \
	materials/Homogeneous.py \
	materials/Material.py \
	materials/MaxwellIsotropic3D.py \
	meshio/__init__.py \
	meshio/CellFilter.py \
	meshio/CellFilterAvg.py \
	meshio/DataWriter.py \
	meshio/DataWriterVTK.py \
	meshio/MeshIOObj.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOCubit.py \
	meshio/MeshIOLagrit.py \
	meshio/OutputDirichlet.py \
	meshio/OutputManager.py \
	meshio/OutputSoln.py \
	meshio/OutputFaultKin.py \
	meshio/OutputMatElastic.py \
	meshio/OutputNeumann.py \
	meshio/OutputSolnSubset.py \
	meshio/SingleOutput.py \
	meshio/VertexFilter.py \
	meshio/VertexFilterVecNorm.py \
	mpi/__init__.py \
	mpi/Communicator.py \
	problems/__init__.py \
	problems/Explicit.py \
	problems/Formulation.py \
	problems/Implicit.py \
	problems/Problem.py \
	problems/TimeDependent.py \
	problems/TimeStep.py \
	problems/TimeStepAdapt.py \
	problems/TimeStepUniform.py \
	problems/TimeStepUser.py \
	solver/__init__.py \
	solver/Solver.py \
	solver/SolverLinear.py \
	solver/SolverNonlinear.py \
	topology/__init__.py \
	topology/Distributor.py \
	topology/Mesh.py \
	topology/SubMesh.py \
	topology/Field.py \
	topology/Fields.py \
	topology/SolutionFields.py \
	topology/Jacobian.py \
	topology/MeshGenerator.py \
	topology/MeshImporter.py \
	topology/MeshRefiner.py \
	topology/RefineUniform.py \
	utils/__init__.py \
	utils/CheckpointTimer.py \
	utils/CppData.py \
	utils/EmptyBin.py \
	utils/EventLogger.py \
	utils/PetscManager.py \
	utils/VTKDataReader.py \
	utils/importing.py \
	utils/profiling.py \
	utils/testarray.py


if ENABLE_TETGEN
  nobase_pkgpyexec_PYTHON += topology/MeshGenSimple.py
endif


# End of file 
