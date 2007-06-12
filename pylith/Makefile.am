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
	bc/BCIntegrator.py \
	bc/BoundaryCondition.py \
	bc/Dirichlet.py \
	bc/BCSingle.py \
	bc/BCTwoSides.py \
	bc/BCFourSides.py \
	bc/BCSixSides.py \
	faults/__init__.py \
	faults/BruneSlipFn.py \
	faults/EqKinSrc.py \
	faults/Fault.py \
	faults/FaultCohesive.py \
	faults/FaultCohesiveKin.py \
	faults/SingleFault.py \
	faults/SlipTimeFn.py \
	feassemble/__init__.py \
	feassemble/Constraint.py \
	feassemble/ElasticityExplicit.py \
	feassemble/ElasticityImplicit.py \
	feassemble/FIATLagrange.py \
	feassemble/FIATSimplex.py \
	feassemble/Integrator.py \
	feassemble/IntegratorElasticity.py \
	feassemble/ReferenceCell.py \
	feassemble/quadrature/__init__.py \
	feassemble/quadrature/Quadrature.py \
	feassemble/quadrature/Quadrature0D.py \
	feassemble/quadrature/Quadrature1D.py \
	feassemble/quadrature/Quadrature1Din2D.py \
	feassemble/quadrature/Quadrature1Din3D.py \
	feassemble/quadrature/Quadrature2D.py \
	feassemble/quadrature/Quadrature2Din3D.py \
	feassemble/quadrature/Quadrature3D.py \
	feassemble/geometry/__init__.py \
	feassemble/geometry/CellGeometry.py \
	feassemble/geometry/GeometryPoint1D.py \
	feassemble/geometry/GeometryPoint2D.py \
	feassemble/geometry/GeometryPoint3D.py \
	feassemble/geometry/GeometryLine1D.py \
	feassemble/geometry/GeometryLine2D.py \
	feassemble/geometry/GeometryLine3D.py \
	feassemble/geometry/GeometryTri2D.py \
	feassemble/geometry/GeometryTri3D.py \
	feassemble/geometry/GeometryQuad2D.py \
	feassemble/geometry/GeometryQuad3D.py \
	feassemble/geometry/GeometryTet3D.py \
	feassemble/geometry/GeometryHex3D.py \
	materials/__init__.py \
	materials/BiMaterial.py \
	materials/ElasticMaterial.py \
	materials/ElasticIsotropic3D.py \
	materials/ElasticPlaneStrain.py \
	materials/ElasticPlaneStress.py \
	materials/ElasticStrain1D.py \
	materials/ElasticStress1D.py \
	materials/Homogeneous.py \
	materials/Material.py \
	materials/MaxwellIsotropic3D.py \
	meshio/__init__.py \
	meshio/MeshIO.py \
	meshio/MeshIOAscii.py \
	meshio/MeshIOCubit.py \
	meshio/MeshIOLagrit.py \
	meshio/SingleOutput.py \
	meshio/SolutionIO.py \
	meshio/SolutionIOVTK.py \
	problems/__init__.py \
	problems/Explicit.py \
	problems/Formulation.py \
	problems/Implicit.py \
	problems/Problem.py \
	problems/TimeDependent.py \
	solver/__init__.py \
	solver/Solver.py \
	solver/SolverLinear.py \
	topology/__init__.py \
	topology/Distributor.py \
	topology/FieldsManager.py \
	topology/Mesh.py \
	topology/MeshGenerator.py \
	topology/MeshImporter.py \
	utils/__init__.py \
	utils/CheckpointTimer.py \
	utils/CppData.py \
	utils/ObjectBin.py \
	utils/PetscManager.py \
	utils/importing.py \
	utils/testarray.py


if ENABLE_TETGEN
  nobase_pkgpyexec_PYTHON += topology/MeshGenSimple.py
endif


# End of file 
