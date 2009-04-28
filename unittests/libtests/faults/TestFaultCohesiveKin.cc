// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveKin.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "data/CohesiveKinData.hh" // USES CohesiveKinData

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKin );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKin::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
  const int nsrcs = 1;
  _eqsrcs.resize(nsrcs);
  _eqsrcs[0] = new EqKinSrc();
  _slipfns.resize(nsrcs);
  _slipfns[0] = new BruneSlipFn();

  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveKin::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  int nsrcs = _eqsrcs.size();
  for (int i=0; i < nsrcs; ++i)
    delete _eqsrcs[i];
  nsrcs = _slipfns.size();
  for (int i=0; i < nsrcs; ++i)
    delete _slipfns[i];
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveKin::testConstructor(void)
{ // testConstructor
  FaultCohesiveKin fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test eqsrc().
void
pylith::faults::TestFaultCohesiveKin::testEqsrc(void)
{ // testEqsrc
  FaultCohesiveKin fault;

  EqKinSrc eqsrcA;
  EqKinSrc eqsrcB;
  const int nsrcs = 2;
  EqKinSrc** sources = new EqKinSrc*[2];
  const char* names[] = {"one", "two"};
  sources[0] = &eqsrcA;
  sources[1] = &eqsrcB;
  fault.eqsrcs(names, sources, nsrcs);
  CPPUNIT_ASSERT(&eqsrcA == fault._eqSrcs["one"]);
  CPPUNIT_ASSERT(&eqsrcB == fault._eqSrcs["two"]);
  delete[] sources; sources = 0;
} // testEqsrc

// ----------------------------------------------------------------------
/// Test needNewJacobian()
void
pylith::faults::TestFaultCohesiveKin::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault.needNewJacobian());
  fault._needNewJacobian = false;
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
} // testNeedNewJacobian

// ----------------------------------------------------------------------
/// Test useSolnIncr()
void
pylith::faults::TestFaultCohesiveKin::testUseSolnIncr(void)
{ // testUseSolnIncr
  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(false, fault._useSolnIncr);
  fault.useSolnIncr(true);
  CPPUNIT_ASSERT_EQUAL(true, fault._useSolnIncr);
} // testUseSolnIncr

// ----------------------------------------------------------------------
// Test useLagrangeConstraints().
void
pylith::faults::TestFaultCohesiveKin::testUseLagrangeConstraints(void)
{ // testUseLagrangeConstraints
  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault._useLagrangeConstraints());
} // testUseLagrangeConstraints

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveKin::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  SieveSubMesh::renumbering_type& renumbering = 
    faultSieveMesh->getRenumbering();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    CPPUNIT_ASSERT(renumbering.find(_data->constraintVertices[iVertex]) !=
		   renumbering.end());
    CPPUNIT_ASSERT_EQUAL(renumbering[_data->constraintVertices[iVertex]],
			 *v_iter);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintVert, iVertex);

  // Check orientation
  const ALE::Obj<RealSection>& orientationSection = 
    fault._fields->get("orientation").section();
  CPPUNIT_ASSERT(!orientationSection.isNull());
  const int cellDim = _data->cellDim;
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = orientationSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(orientationSize, fiberDim);
    const double* orientationVertex =
      orientationSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != orientationVertex);

    const double tolerance = 1.0e-06;
    for (int i=0; i < orientationSize; ++i) {
      const int index = iVertex*orientationSize+i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[index],
				   orientationVertex[i], tolerance);
    } // for
  } // for

  // Check area
  const ALE::Obj<RealSection>& areaSection =
    fault._fields->get("area").section();
  CPPUNIT_ASSERT(!areaSection.isNull());
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = areaSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(1, fiberDim);
    const double* areaVertex = areaSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != areaVertex);

    const double tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->area[iVertex], areaVertex[0],
				 tolerance);
  } // for

  // Check pseudoStiffness
  const ALE::Obj<RealSection>& stiffnessSection =
    fault._fields->get("pseudostiffness").section();
  CPPUNIT_ASSERT(!stiffnessSection.isNull());
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = stiffnessSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(1, fiberDim);
    const double* stiffnessVertex = stiffnessSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != stiffnessVertex);
    const double tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->pseudoStiffness, stiffnessVertex[0],
				 tolerance);
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidual(void)
{ // testIntegrateResidual
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());

  const ALE::Obj<RealSection>& solutionSection = 
    fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  
  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  { // Integrate residual with solution (as opposed to solution increment).
    fault.useSolnIncr(false);
    fault.integrateResidual(residual, t, &fields);

    //residual->view("RESIDUAL"); // DEBUGGING

    // Check values
    const double* valsE = _data->valsResidual;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      const bool isConstraint = _isConstraintVertex(*v_iter);
      if (!isConstraint) {
	const double pseudoStiffness = _data->pseudoStiffness;
	for (int i=0; i < fiberDimE; ++i) {
	  const int index = iVertex*spaceDim+i;
	  const double valE = valsE[index] * pseudoStiffness;
	  if (valE > tolerance)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
	  else
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
	} // for
      } else {
	const double valE = 0.0; // no contribution
	for (int i=0; i < fiberDimE; ++i)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Integrate residual with solution (as opposed to solution increment).

  residual.zero();
  { // Integrate residual with solution increment.
    fault.useSolnIncr(true);
    fault.integrateResidual(residual, t, &fields);

    residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const double* valsE = _data->valsResidualIncr;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      const bool isConstraint = _isConstraintVertex(*v_iter);
      if (!isConstraint) {
	const double pseudoStiffness = _data->pseudoStiffness;
	for (int i=0; i < fiberDimE; ++i) {
	  const int index = iVertex*spaceDim+i;
	  const double valE = valsE[index] * pseudoStiffness;
	  if (valE > tolerance)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
	  else
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
	} // for
      } else {
	const double valE = 0.0; // no contribution
	for (int i=0; i < fiberDimE; ++i)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Integrate residual with solution increment.
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const ALE::Obj<RealSection>& solutionSection = fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  const int spaceDim = _data->spaceDim;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  } // for
  
  topology::Jacobian jacobian(fields);

  const double t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  jacobian.assemble("final_assembly");

  //MatView(jacobian, PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  const double* valsE = _data->valsJacobian;
  const int nrowsE = solutionSection->sizeWithBC();
  const int ncolsE = nrowsE;

  int nrows = 0;
  int ncols = 0;
  PetscMat jacobianMat = jacobian.matrix();
  MatGetSize(jacobianMat, &nrows, &ncols);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  PetscMat jSparseAIJ;
  MatConvert(jacobianMat, MATSEQAIJ, MAT_INITIAL_MATRIX, &jSparseAIJ);
  MatConvert(jSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);

  double_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const double tolerance = 1.0e-06;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      const double valE = 0.0;
#if 0 // DEBUGGING
      if (fabs(valE-vals[index]) > tolerance)
	std::cout << "ERROR: iRow: " << iRow << ", iCol: " << iCol
		  << "valE: " << valE
		  << ", val: " << vals[index]
		  << std::endl;
#endif // DEBUGGING
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[index], tolerance);
    } // for
  MatDestroy(jDense);
  MatDestroy(jSparseAIJ);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateResidualAssembled().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidualAssembled(void)
{ // testIntegrateResidualAssembled
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());

  const ALE::Obj<RealSection>& solutionSection = fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  
  const double t = 2.134;
  const double dt = 0.01;
  fault.timeStep(dt);
  { // Integrate residual with solution (as opposed to solution increment).
    fault.useSolnIncr(false);
    fault.integrateResidualAssembled(residual, t, &fields);

    //residual->view("RESIDUAL"); // DEBUGGING

    // Check values
    const double* valsE = _data->valsResidual;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      const bool isConstraint = _isConstraintVertex(*v_iter);
      if (!isConstraint) {
	const double valE = 0.0;
	for (int i=0; i < fiberDimE; ++i)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } else {
	const double pseudoStiffness = _data->pseudoStiffness;
	for (int i=0; i < fiberDimE; ++i) {
	  const int index = iVertex*spaceDim+i;
	  const double valE = valsE[index] * pseudoStiffness;
	  if (valE > tolerance)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
	  else
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
	} // for
      } // if/else
    } // for
  } // Integrate residual with solution (as opposed to solution increment).

  residual.zero();
  { // Integrate residual with solution increment.
    fault.useSolnIncr(true);
    fault.integrateResidualAssembled(residual, t, &fields);

    //residual->view("RESIDUAL"); // DEBUGGING

    // Check values
    const double* valsE = _data->valsResidualIncr;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const double tolerance = 1.0e-06;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const double* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      const bool isConstraint = _isConstraintVertex(*v_iter);
      if (!isConstraint) {
	const double valE = 0.0;
	for (int i=0; i < fiberDimE; ++i)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } else {
	const double pseudoStiffness = _data->pseudoStiffness;
	for (int i=0; i < fiberDimE; ++i) {
	  const int index = iVertex*spaceDim+i;
	  const double valE = valsE[index] * pseudoStiffness;
	  if (valE > tolerance)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
	  else
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
	} // for
      } // if/else
    } // for
  } // Integrate residual with solution increment.
} // testIntegrateResidualAssembled

// ----------------------------------------------------------------------
// Test integrateJacobianAssembled().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobianAssembled(void)
{ // testIntegrateJacobianAssembled
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  const ALE::Obj<RealSection>& solutionSection = fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex)
    solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);

  topology::Jacobian jacobian(fields);

  const double t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  jacobian.assemble("final_assembly");

  //MatView(jacobian, PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  const double* valsE = _data->valsJacobian;
  const int nrowsE = solutionSection->sizeWithBC();
  const int ncolsE = nrowsE;

  PetscMat jacobianMat = jacobian.matrix();

  int nrows = 0;
  int ncols = 0;
  MatGetSize(jacobianMat, &nrows, &ncols);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  PetscMat jSparseAIJ;
  MatConvert(jacobianMat, MATSEQAIJ, MAT_INITIAL_MATRIX, &jSparseAIJ);
  MatConvert(jSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);

  double_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const double tolerance = 1.0e-06;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      const double pseudoStiffness = 
	(iRow > iCol) ? 1.0 : _data->pseudoStiffness;
      const double valE = valsE[index] * pseudoStiffness;
#if 0 // DEBUGGING
      if (fabs(valE-vals[index]) > tolerance)
	std::cout << "ERROR: iRow: " << iRow << ", iCol: " << iCol
		  << "valE: " << valE
		  << ", val: " << vals[index]
		  << std::endl;
#endif // DEBUGGING
      if (fabs(valE) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valE, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[index], tolerance);
    } // for
  MatDestroy(jDense);
  MatDestroy(jSparseAIJ);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
} // testIntegrateJacobianAssembled

// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::faults::TestFaultCohesiveKin::testUpdateStateVars(void)
{ // testUpdateStateVars
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  const ALE::Obj<RealSection>& solutionSection = fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());
  { // setup solution
    solutionSection->zero();
    
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
    int iVertex = 0;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex)
      solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  } // setup solution
  topology::Field<topology::Mesh>& residual = fields.get("residual");

  const double t = 2.134;
  const double dt = 0.01;
  fault.useSolnIncr(false);
  fault.timeStep(dt);
  fault.integrateResidualAssembled(residual, t, &fields);
  fault.updateStateVars(t, &fields);

  CPPUNIT_ASSERT(0 != fault._faultMesh);
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();

  // Compute expected cumulative slip using eqsrcs
  topology::Field<topology::SubMesh> cumSlipE(*fault._faultMesh);
  cumSlipE.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  cumSlipE.allocate();
  const ALE::Obj<RealSection> cumSlipESection = cumSlipE.section();
  CPPUNIT_ASSERT(!cumSlipESection.isNull());

  const ALE::Obj<RealSection> cumSlipSection =
    fault._fields->get("cumulative slip").section();
  CPPUNIT_ASSERT(!cumSlipSection.isNull());

  const FaultCohesiveKin::srcs_type::const_iterator srcsEnd = fault._eqSrcs.end();
  for (FaultCohesiveKin::srcs_type::iterator s_iter=fault._eqSrcs.begin(); 
       s_iter != srcsEnd; 
       ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    if (t >= src->originTime())
      src->slip(&cumSlipE, t);
  } // for

  int iVertex = 0;
  const double tolerance = 1.0e-06;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const SieveSubMesh::point_type meshVertex = _data->constraintVertices[iVertex];
    bool found = false;
    for(SieveSubMesh::renumbering_type::const_iterator r_iter = renumbering.begin();
	r_iter != renumbering.end();
	++r_iter) {
      if (r_iter->second == *v_iter) {
        found = true;
        break;
      } // if
    } // for
    CPPUNIT_ASSERT(found);

    // Check _cumSlip
    int fiberDim = cumSlipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* slipV = cumSlipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipV);

    const double* slipE = cumSlipESection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipE);

    for (int iDim=0; iDim < spaceDim; ++iDim) {
      if (slipE[iDim] > 1.0) 
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, slipV[iDim]/slipE[iDim], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE[iDim], slipV[iDim], tolerance);
    } // for
  } // for
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcTractionsChange().
void
pylith::faults::TestFaultCohesiveKin::testCalcTractionsChange(void)
{ // testCalcTractionsChange
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  
  const int spaceDim = _data->spaceDim;
  const ALE::Obj<RealSection>& solutionSection = fields.get("solution").section();
  CPPUNIT_ASSERT(!solutionSection.isNull());
  { // setup solution
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
    int iVertex = 0;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      solutionSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
    } // for
  } // setup solution

  CPPUNIT_ASSERT(0 != fault._faultMesh);
  topology::Field<topology::SubMesh> tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();
  const ALE::Obj<RealSection>& tractionsSection = tractions.section();
  CPPUNIT_ASSERT(!tractionsSection.isNull());

  const double t = 0;
  fault.updateStateVars(t, &fields);  
  fault._calcTractionsChange(&tractions, fields.get("solution"));

  int iVertex = 0;
  const double tolerance = 1.0e-06;
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering = 
    faultSieveMesh->getRenumbering();
  const SieveMesh::renumbering_type::const_iterator renumberingBegin = 
    renumbering.begin();
  const SieveMesh::renumbering_type::const_iterator renumberingEnd = 
    renumbering.end();
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    SieveMesh::point_type meshVertex = -1;
    bool found = false;

    for (SieveMesh::renumbering_type::const_iterator r_iter = renumberingBegin;
	 r_iter != renumberingEnd;
	 ++r_iter) {
      if (r_iter->second == *v_iter) {
        meshVertex = r_iter->first;
        found = true;
        break;
      } // if
    } // for
    CPPUNIT_ASSERT(found);
    int fiberDim = tractionsSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* tractionsVertex = tractionsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != tractionsVertex);

    fiberDim = solutionSection->getFiberDimension(meshVertex);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const double* solutionVertex = solutionSection->restrictPoint(meshVertex);
    CPPUNIT_ASSERT(0 != solutionVertex);

    const double scale = _data->pseudoStiffness / _data->area[iVertex];
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const double tractionE = solutionVertex[iDim] * scale;
      if (tractionE > 1.0) 
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tractionsVertex[iDim]/tractionE,
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, tractionsVertex[iDim],
				     tolerance);
    } // for
  } // for
} // testCalcTractionsChange

// ----------------------------------------------------------------------
// Initialize FaultCohesiveKin interface condition.
void
pylith::faults::TestFaultCohesiveKin::_initialize(
					topology::Mesh* const mesh,
					FaultCohesiveKin* const fault,
					topology::SolutionFields* const fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  
  //(*mesh)->setDebug(true); // DEBUGGING
  
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup earthquake source
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
  ioFinalSlip.filename(_data->finalSlipFilename);
  dbFinalSlip.ioHandler(&ioFinalSlip);
  
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
  ioSlipTime.filename(_data->slipTimeFilename);
  dbSlipTime.ioHandler(&ioSlipTime);
  
  spatialdata::spatialdb::SimpleDB dbRiseTime("rise time");
  spatialdata::spatialdb::SimpleIOAscii ioRiseTime;
  ioRiseTime.filename(_data->riseTimeFilename);
  dbRiseTime.ioHandler(&ioRiseTime);
  
  const int nsrcs = _eqsrcs.size();
  CPPUNIT_ASSERT(nsrcs == _slipfns.size());
  EqKinSrc** sources = new EqKinSrc*[nsrcs];
  char** names = new char*[nsrcs];
  for (int i=0; i < nsrcs; ++i) {
    _slipfns[i]->dbFinalSlip(&dbFinalSlip);
    _slipfns[i]->dbSlipTime(&dbSlipTime);
    _slipfns[i]->dbRiseTime(&dbRiseTime);
    
    _eqsrcs[i]->slipfn(_slipfns[i]);
    sources[i] = _eqsrcs[i];
    names[i] = new char[2];
    names[i][0] = 'a' + i;
    names[i][1] = '\0';
  } // for
  
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->eqsrcs(const_cast<const char**>(names), sources, nsrcs);
  fault->adjustTopology(mesh, _flipFault);
  
  const double upDir[] = { 0.0, 0.0, 1.0 };
  const double normalDir[] = { 1.0, 0.0, 0.0 };
  
  spatialdata::spatialdb::SimpleDB dbMatProp("material properties");
  spatialdata::spatialdb::SimpleIOAscii ioMatProp;
  ioMatProp.filename(_data->matPropsFilename);
  dbMatProp.ioHandler(&ioMatProp);
  
  fault->initialize(*mesh, upDir, normalDir, &dbMatProp); 
  
  delete[] sources; sources = 0;
  for (int i=0; i < nsrcs; ++i)
    delete[] names[i];
  delete[] names; names = 0;
  
  // Setup fields
  fields->add("residual", "residual");
  fields->add("solution", "displacement");
  fields->solutionName("solution");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  fields->copyLayout("residual");
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveKin::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  assert(0 != _data);

  const int numConstraintVert = _data->numConstraintVert;
  bool isFound = false;
  for (int i=0; i < _data->numConstraintVert; ++i)
    if (_data->constraintVertices[i] == vertex) {
      isFound = true;
      break;
    } // if
  return isFound;
} // _isConstraintVertex


// End of file 
