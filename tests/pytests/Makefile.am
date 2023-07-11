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

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	apps/__init__.py \
	apps/TestPetscApplication.py \
	apps/TestPyLithApp.py \
	apps/TestEqInfoApp.py \
	bc/__init__.py \
	bc/TestDirichletTimeDependent.py \
	bc/TestNeumannTimeDependent.py \
	bc/TestAbsorbingDampers.py \
	bc/TestAuxiliarySubfields.py \
	bc/TestZeroDB.py \
	faults/__init__.py \
	faults/TestFaultCohesive.py \
	faults/TestFaultCohesiveKin.py \
	faults/TestFaultCohesiveImpulses.py \
	faults/TestKinSrc.py \
	faults/TestKinSrcConstRate.py \
	faults/TestKinSrcStep.py \
	faults/TestKinSrcRamp.py \
	faults/TestKinSrcBrune.py \
	faults/TestKinSrcLiuCos.py \
	faults/TestKinSrcTimeHistory.py \
	faults/TestSingleRupture.py \
	materials/__init__.py \
	materials/TestAuxiliarySubfields.py \
	materials/TestDerivedSubfields.py \
	materials/TestElasticity.py \
	materials/TestHomogeneous.py \
	materials/TestIncompressibleElasticity.py \
	materials/TestMaterial.py \
	materials/TestPoroelasticity.py \
	materials/TestRheologies.py \
	materials/TestRheologyLinearElasticity.py \
	materials/TestRheologyPoroelasticity.py \
	materials/TestRheologyViscoelasticity.py \
	meshio/__init__.py \
	meshio/TestDataWriter.py \
	meshio/TestDataWriterHDF5.py \
	meshio/TestDataWriterHDF5Ext.py \
	meshio/TestDataWriterVTK.py \
	meshio/TestMeshIOAscii.py \
	meshio/TestMeshIOCubit.py \
	meshio/TestOutputManagerMesh.py \
	meshio/TestOutputManagerSubmesh.py \
	meshio/TestOutputObserver.py \
	meshio/TestOutputPhysics.py \
	meshio/TestOutputSoln.py \
	meshio/TestOutputSolnBoundary.py \
	meshio/TestOutputSolnDomain.py \
	meshio/TestOutputSolnPoints.py \
	meshio/TestOutputSolnSubset.py \
	meshio/TestOutputTrigger.py \
	meshio/TestOutputTriggerStep.py \
	meshio/TestOutputTriggerTime.py \
	meshio/TestPointsList.py \
	meshio/TestSingleOutput.py \
	meshio/TestXdmf.py \
	mpi/__init__.py \
	mpi/TestCommunicator.py \
	mpi/TestReduce.py \
	problems/__init__.py \
	problems/TestInitialCondition.py \
	problems/TestInitialConditionDomain.py \
	problems/TestInitialConditionPatch.py \
	problems/TestPhysics.py \
	problems/TestProblem.py \
	problems/TestProblemDefaults.py \
	problems/TestProgressMonitor.py \
	problems/TestProgressMonitorTime.py \
	problems/TestSingleObserver.py \
	problems/TestSolnDisp.py \
	problems/TestSolnDispLagrange.py \
	problems/TestSolnDispPres.py \
	problems/TestSolnDispPresLagrange.py \
	problems/TestSolnDispTracStrain.py \
	problems/TestSolnDispVel.py \
	problems/TestSolnDispVelLagrange.py \
	problems/TestSolution.py \
	problems/TestSolutionSubfield.py \
	problems/TestSolutionSubfields.py \
	problems/TestSubfieldDisplacement.py \
	problems/TestSubfieldLagrangeFault.py \
	problems/TestSubfieldPressure.py \
	problems/TestSubfieldTemperature.py \
	problems/TestSubfieldTraceStrain.py \
	problems/TestSubfieldVelocity.py \
	problems/TestTimeDependent.py \
	topology/__init__.py \
	topology/TestDistributor.py \
	topology/TestField.py \
	topology/TestMesh.py \
	topology/TestMeshGenerator.py \
	topology/TestMeshImporter.py \
	topology/TestMeshRefiner.py \
	topology/TestRefineUniform.py \
	topology/TestReverseCuthillMcKee.py \
	topology/TestSubfield.py \
	utils/__init__.py \
	utils/TestCollectVersionInfo.py \
	utils/TestConstants.py \
	utils/TestDependenciesVersion.py \
	utils/TestDumpParameters.py \
	utils/TestDumpParametersAscii.py \
	utils/TestDumpParametersJson.py \
	utils/TestEmptyBin.py \
	utils/TestEventLogger.py \
	utils/TestNullComponent.py \
	utils/TestPetscManager.py \
	utils/TestPetscVersion.py \
	utils/TestProfiling.py \
	utils/TestPylithVersion.py

dist_noinst_DATA = \
	meshio/data/tri3.h5 \
	meshio/data/tri3.xmf \
	meshio/data/tri3_vertex.h5 \
	meshio/data/tri3_vertex.xmf \
	meshio/data/tri3_cell.h5 \
	meshio/data/tri3_cell.xmf \
	meshio/data/tri3_points.h5 \
	meshio/data/tri3_points.xmf \
	meshio/data/tri3_points_vertex.h5 \
	meshio/data/tri3_points_vertex.xmf \
	meshio/data/tri3_surf.h5 \
	meshio/data/tri3_surf.xmf \
	meshio/data/tri3_surf_vertex.h5 \
	meshio/data/tri3_surf_vertex.xmf \
	meshio/data/tri3_surf_cell.h5 \
	meshio/data/tri3_surf_cell.xmf \
	meshio/data/quad4.h5 \
	meshio/data/quad4.xmf \
	meshio/data/quad4_vertex.h5 \
	meshio/data/quad4_vertex.xmf \
	meshio/data/quad4_cell.h5 \
	meshio/data/quad4_cell.xmf \
	meshio/data/quad4_points.h5 \
	meshio/data/quad4_points.xmf \
	meshio/data/quad4_points_vertex.h5 \
	meshio/data/quad4_points_vertex.xmf \
	meshio/data/quad4_surf.h5 \
	meshio/data/quad4_surf.xmf \
	meshio/data/quad4_surf_vertex.h5 \
	meshio/data/quad4_surf_vertex.xmf \
	meshio/data/quad4_surf_cell.h5 \
	meshio/data/quad4_surf_cell.xmf \
	meshio/data/tet4.h5 \
	meshio/data/tet4.xmf \
	meshio/data/tet4_vertex.h5 \
	meshio/data/tet4_vertex.xmf \
	meshio/data/tet4_cell.h5 \
	meshio/data/tet4_cell.xmf \
	meshio/data/tet4_points.h5 \
	meshio/data/tet4_points.xmf \
	meshio/data/tet4_points_vertex.h5 \
	meshio/data/tet4_points_vertex.xmf \
	meshio/data/tet4_surf.h5 \
	meshio/data/tet4_surf.xmf \
	meshio/data/tet4_surf_vertex.h5 \
	meshio/data/tet4_surf_vertex.xmf \
	meshio/data/tet4_surf_cell.h5 \
	meshio/data/tet4_surf_cell.xmf \
	meshio/data/hex8.h5 \
	meshio/data/hex8.xmf \
	meshio/data/hex8_vertex.h5 \
	meshio/data/hex8_vertex.xmf \
	meshio/data/hex8_cell.h5 \
	meshio/data/hex8_cell.xmf \
	meshio/data/hex8_points.h5 \
	meshio/data/hex8_points.xmf \
	meshio/data/hex8_points_vertex.h5 \
	meshio/data/hex8_points_vertex.xmf \
	meshio/data/hex8_surf.h5 \
	meshio/data/hex8_surf.xmf \
	meshio/data/hex8_surf_vertex.h5 \
	meshio/data/hex8_surf_vertex.xmf \
	meshio/data/hex8_surf_cell.h5 \
	meshio/data/hex8_surf_cell.xmf

noinst_TMP = \
	meshio/tri3.xmf \
	meshio/tri3_vertex.xmf \
	meshio/tri3_cell.xmf \
	meshio/tri3_points.xmf \
	meshio/tri3_points_vertex.xmf \
	meshio/tri3_surf.xmf \
	meshio/tri3_surf_vertex.xmf \
	meshio/tri3_surf_cell.xmf \
	meshio/quad4.xmf \
	meshio/quad4_vertex.xmf \
	meshio/quad4_cell.xmf \
	meshio/quad4_points.xmf \
	meshio/quad4_points_vertex.xmf \
	meshio/quad4_surf.xmf \
	meshio/quad4_surf_vertex.xmf \
	meshio/quad4_surf_cell.xmf \
	meshio/tet4.xmf \
	meshio/tet4_vertex.xmf \
	meshio/tet4_cell.xmf \
	meshio/tet4_points.xmf \
	meshio/tet4_points_vertex.xmf \
	meshio/tet4_surf.xmf \
	meshio/tet4_surf_vertex.xmf \
	meshio/tet4_surf_cell.xmf \
	meshio/hex8.xmf \
	meshio/hex8_vertex.xmf \
	meshio/hex8_cell.xmf \
	meshio/hex8_points.xmf \
	meshio/hex8_points_vertex.xmf \
	meshio/hex8_surf.xmf \
	meshio/hex8_surf_vertex.xmf \
	meshio/hex8_surf_cell.xmf


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file 
