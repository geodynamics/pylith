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
# Copyright (c) 2010-2011 University of California, Davis
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
	apps/PrepMeshApp.py \
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
	faults/Fault.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveKin.py \
	faults/FaultCohesiveDyn.py \
	faults/FaultCohesiveTract.py \
	faults/LiuCosSlipFn.py \
	faults/SingleRupture.py \
	faults/SlipTimeFn.py \
	faults/StepSlipFn.py \
	faults/TimeHistorySlipFn.py \
	feassemble/__init__.py \
	feassemble/Constraint.py \
	feassemble/ElasticityExplicit.py \
	feassemble/ElasticityExplicitTet4.py \
	feassemble/ElasticityExplicitTri3.py \
	feassemble/ElasticityExplicitLgDeform.py \
	feassemble/ElasticityImplicit.py \
	feassemble/ElasticityImplicitLgDeform.py \
	feassemble/FIATLagrange.py \
	feassemble/FIATSimplex.py \
	feassemble/Integrator.py \
	feassemble/IntegratorElasticity.py \
	feassemble/IntegratorElasticityLgDeform.py \
	feassemble/ReferenceCell.py \
	feassemble/Quadrature.py \
	feassemble/CellGeometry.py \
	friction/__init__.py \
	friction/FrictionModel.py \
	friction/StaticFriction.py \
	friction/SlipWeakening.py \
	friction/RateStateAgeing.py \
	friction/TimeWeakening.py \
	materials/__init__.py \
	materials/ElasticMaterial.py \
	materials/ElasticIsotropic3D.py \
	materials/ElasticPlaneStrain.py \
	materials/ElasticPlaneStress.py \
	materials/ElasticStrain1D.py \
	materials/ElasticStress1D.py \
	materials/GenMaxwellIsotropic3D.py \
	materials/GenMaxwellPlaneStrain.py \
	materials/GenMaxwellQpQsIsotropic3D.py \
	materials/Homogeneous.py \
	materials/Material.py \
	materials/MaxwellIsotropic3D.py \
	materials/MaxwellPlaneStrain.py \
	materials/PowerLaw3D.py \
	materials/DruckerPrager3D.py \
	meshio/__init__.py \
	meshio/CellFilter.py \
	meshio/CellFilterAvgMesh.py \
	meshio/CellFilterAvgSubMesh.py \
	meshio/DataWriter.py \
	meshio/DataWriterVTK.py \
	meshio/DataWriterVTKMesh.py \
	meshio/DataWriterVTKSubMesh.py \
	meshio/DataWriterVTKSubSubMesh.py \
	meshio/MeshIOObj.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOLagrit.py \
	meshio/MeshIOSieve.py \
	meshio/OutputDirichlet.py \
	meshio/OutputManager.py \
	meshio/OutputManagerMesh.py \
	meshio/OutputManagerSubMesh.py \
	meshio/OutputSoln.py \
	meshio/OutputFaultKin.py \
	meshio/OutputFaultDyn.py \
	meshio/OutputMatElastic.py \
	meshio/OutputNeumann.py \
	meshio/OutputSolnSubset.py \
	meshio/SingleOutput.py \
	meshio/VertexFilter.py \
	meshio/VertexFilterVecNormMesh.py \
	meshio/VertexFilterVecNormSubMesh.py \
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
	problems/ExplicitLumped.py \
	problems/ExplicitLumpedTri3.py \
	problems/ExplicitLumpedTet4.py \
	problems/ExplicitLgDeform.py \
	problems/Formulation.py \
	problems/Implicit.py \
	problems/ImplicitLgDeform.py \
	problems/Problem.py \
	problems/Solver.py \
	problems/SolverLinear.py \
	problems/SolverNonlinear.py \
	problems/SolverLumped.py \
	problems/TimeDependent.py \
	problems/TimeStep.py \
	problems/TimeStepAdapt.py \
	problems/TimeStepUniform.py \
	problems/TimeStepUser.py \
	topology/__init__.py \
	topology/Distributor.py \
	topology/Mesh.py \
	topology/SubMesh.py \
	topology/Field.py \
	topology/Fields.py \
	topology/SolutionFields.py \
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
	tests/PhysicalProperties.py \
	tests/Solution.py \
	tests/StateVariables.py \
	tests/Fault.py

if ENABLE_CUBIT
  nobase_pkgpyexec_PYTHON += meshio/MeshIOCubit.py
endif

if ENABLE_HDF5
  nobase_pkgpyexec_PYTHON += \
	meshio/DataWriterHDF5.py \
	meshio/DataWriterHDF5Mesh.py \
	meshio/DataWriterHDF5SubMesh.py \
	meshio/DataWriterHDF5SubSubMesh.py \
	meshio/DataWriterHDF5Ext.py \
	meshio/DataWriterHDF5ExtMesh.py \
	meshio/DataWriterHDF5ExtSubMesh.py \
	meshio/DataWriterHDF5ExtSubSubMesh.py \
	meshio/Xdmf.py
endif

if ENABLE_CUDA
  nobase_pkgpyexec_PYTHON += \
	feassemble/ElasticityImplicitCUDA.py \
	problems/ImplicitCUDA.py
endif

if ENABLE_TETGEN
  nobase_pkgpyexec_PYTHON += topology/MeshGenSimple.py
endif


# End of file 
