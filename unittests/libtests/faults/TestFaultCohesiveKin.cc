// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
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
  fault.eqsrcs(names, nsrcs, sources, nsrcs);
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
  CPPUNIT_ASSERT_EQUAL(true, fault.useLagrangeConstraints());
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
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  int iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    CPPUNIT_ASSERT(renumbering.find(_data->verticesLagrange[iVertex]) !=
		   renumbering.end());
#if 0
    CPPUNIT_ASSERT_EQUAL(renumbering[_data->verticesLagrange[iVertex]],
			 *v_iter);
#endif
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, iVertex);

  // Check cohesive vertex info
  const int numFaultVertices = _data->numFaultVertices;
  CPPUNIT_ASSERT_EQUAL(numFaultVertices, int(fault._cohesiveVertices.size()));
  for (int i=0; i < numFaultVertices; ++i) {
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[i], 
			 fault._cohesiveVertices[i].fault);
    CPPUNIT_ASSERT_EQUAL(_data->verticesLagrange[i], 
			 fault._cohesiveVertices[i].lagrange);
    CPPUNIT_ASSERT_EQUAL(_data->verticesNegative[i], 
			 fault._cohesiveVertices[i].negative);
    CPPUNIT_ASSERT_EQUAL(_data->verticesPositive[i], 
			 fault._cohesiveVertices[i].positive);
  } // for

  // Check cohesive cell info
  const int numCohesiveCells = _data->numCohesiveCells;
  CPPUNIT_ASSERT_EQUAL(numCohesiveCells, int(fault._cohesiveToFault.size()));
  std::map<topology::Mesh::SieveMesh::point_type,
    topology::Mesh::SieveMesh::point_type>::iterator m_iterator = 
    fault._cohesiveToFault.begin();
  std::map<topology::Mesh::SieveMesh::point_type,
    topology::Mesh::SieveMesh::point_type>::const_iterator mapEnd = 
    fault._cohesiveToFault.end();
  for (int i=0; i < numCohesiveCells; ++i, ++m_iterator) {
    CPPUNIT_ASSERT(mapEnd != m_iterator);
    CPPUNIT_ASSERT_EQUAL(_data->cellMappingFault[i], m_iterator->second);
    CPPUNIT_ASSERT_EQUAL(_data->cellMappingCohesive[i], m_iterator->first);
  } // for

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING
  const ALE::Obj<RealSection>& orientationSection = 
    fault._fields->get("orientation").section();
  CPPUNIT_ASSERT(!orientationSection.isNull());
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  iVertex = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = orientationSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(orientationSize, fiberDim);
    const PylithScalar* orientationVertex =
      orientationSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != orientationVertex);

    const PylithScalar tolerance = 1.0e-06;
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
    const PylithScalar* areaVertex = areaSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != areaVertex);

    const PylithScalar tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->area[iVertex], areaVertex[0],
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

  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());

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
    dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  
  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);

  { // Integrate residual with disp (as opposed to disp increment).
    fault.useSolnIncr(false);
    fault.integrateResidual(residual, t, &fields);

    //residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const PylithScalar* valsE = _data->residual;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const PylithScalar* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const PylithScalar valE = valsE[index];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Integrate residual with disp (as opposed to disp increment).

  residual.zero();
  { // Integrate residual with disp increment.
    fault.useSolnIncr(true);
    fault.integrateResidual(residual, t, &fields);

    //residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const PylithScalar* valsE = _data->residualIncr;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
	 v_iter != verticesEnd;
	 ++v_iter, ++iVertex) {
      const int fiberDim = residualSection->getFiberDimension(*v_iter);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
      const PylithScalar* vals = residualSection->restrictPoint(*v_iter);
      CPPUNIT_ASSERT(0 != vals);
      
      for (int i = 0; i < fiberDimE; ++i) {
        const int index = iVertex * spaceDim + i;
        const PylithScalar valE = valsE[index];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[i], tolerance);
      } // for
    } // for
  } // Integrate residual with disp increment.
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

  const int spaceDim = _data->spaceDim;
  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());

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
    dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);

  topology::Jacobian jacobian(fields.solution());

  const PylithScalar t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  jacobian.assemble("final_assembly");

  //MatView(jacobian.matrix(), PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  const PylithScalar* valsE = _data->jacobian;
  const int nrowsE = dispSection->sizeWithBC();
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

  scalar_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const PylithScalar tolerance = 1.0e-06;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      const PylithScalar valE = valsE[index];
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
  MatDestroy(&jDense);
  MatDestroy(&jSparseAIJ);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobian() with lumped Jacobian.
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const PylithScalar t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
  jacobian.complete();

  // jacobian.view("JACOBIAN"); // DEBUGGING

  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  CPPUNIT_ASSERT(!jacobianSection.isNull());

  int iVertex = 0;
  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _data->spaceDim;
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
    int fiberDim = jacobianSection->getFiberDimension(meshVertex);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* jacobianVertex = jacobianSection->restrictPoint(meshVertex);
    CPPUNIT_ASSERT(0 != jacobianVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobianVertex[iDim],
				   tolerance);
  } // for
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
// Test adjustSolnLumped().
void
pylith::faults::TestFaultCohesiveKin::testAdjustSolnLumped(void)
{ // testAdjustSolnLumped
  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  { // setup disp
    const ALE::Obj<RealSection>& dispTSection = fields.get("disp(t)").section();
    CPPUNIT_ASSERT(!dispTSection.isNull());
    int iVertex = 0;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
        dispTSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
    } // for
  } // setup disp

  // compute residual so that slip and residual are setup
  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);
  residual.complete();

  { // setup disp increment
    const ALE::Obj<RealSection>& dispIncrSection = fields.get("dispIncr(t->t+dt)").section();
    CPPUNIT_ASSERT(!dispIncrSection.isNull());
    int iVertex = 0;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
        dispIncrSection->updatePoint(*v_iter, &_data->fieldIncr[iVertex*spaceDim]);
    } // for
  } // setup disp increment

  // Set Jacobian values
  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();
  { // setup disp
    const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
    CPPUNIT_ASSERT(!jacobianSection.isNull());
    int iVertex = 0;
    for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
        v_iter != verticesEnd;
        ++v_iter, ++iVertex) {
        jacobianSection->updatePoint(*v_iter, &_data->jacobianLumped[iVertex*spaceDim]);
    } // for
  } // setup disp
  jacobian.complete();

  topology::Field<topology::Mesh>& solution = fields.get("dispIncr(t->t+dt)");
  fault.adjustSolnLumped(&fields, jacobian);
  const topology::Field<topology::Mesh>& dispIncrAdj = 
    fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  //solution.view("SOLUTION AFTER ADJUSTMENT"); // DEBUGGING

  const ALE::Obj<RealSection>& solutionSection = solution.section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  int i = 0;
  const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  const PylithScalar* solutionE = _data->fieldIncrAdjusted;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    const int fiberDim = solutionSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* solutionVertex = solutionSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != solutionVertex);
    for (int iDim=0; iDim < spaceDim; ++iDim, ++i)
      if (0.0 != solutionE[i])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, solutionVertex[iDim]/solutionE[i],
          tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(solutionE[i], solutionVertex[iDim],
				     tolerance);
  } // for
} // testAdjustSolnLumped

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
  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());
  { // setup disp
    dispSection->zero();
    
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
      dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  } // setup disp
  topology::Field<topology::Mesh>& residual = fields.get("residual");

  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.useSolnIncr(false);
  fault.timeStep(dt);
  fault.integrateResidual(residual, t, &fields);
  fault.updateStateVars(t, &fields);

  CPPUNIT_ASSERT(0 != fault._faultMesh);
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault._faultMesh->sieveMesh();
  CPPUNIT_ASSERT(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();

  // Compute expected slip using eqsrcs
  topology::Field<topology::SubMesh> slipE(*fault._faultMesh);
  slipE.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slipE.allocate();
  const ALE::Obj<RealSection> slipESection = slipE.section();
  CPPUNIT_ASSERT(!slipESection.isNull());

  const ALE::Obj<RealSection> slipSection =
    fault._fields->get("slip").section();
  CPPUNIT_ASSERT(!slipSection.isNull());

  const FaultCohesiveKin::srcs_type::const_iterator srcsEnd = fault._eqSrcs.end();
  for (FaultCohesiveKin::srcs_type::iterator s_iter=fault._eqSrcs.begin(); 
       s_iter != srcsEnd; 
       ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    if (t >= src->originTime())
      src->slip(&slipE, t);
  } // for

  int iVertex = 0;
  const PylithScalar tolerance = 1.0e-06;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++iVertex) {
    const SieveSubMesh::point_type meshVertex = _data->verticesLagrange[iVertex];
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

    // Check _slip
    int fiberDim = slipSection->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* slipV = slipSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != slipV);

    const PylithScalar* slipE = slipESection->restrictPoint(*v_iter);
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
  const ALE::Obj<RealSection>& dispSection = fields.get("disp(t)").section();
  CPPUNIT_ASSERT(!dispSection.isNull());
  { // setup disp
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
      dispSection->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
    } // for
  } // setup disp

  CPPUNIT_ASSERT(0 != fault._faultMesh);
  topology::Field<topology::SubMesh> tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();
  const ALE::Obj<RealSection>& tractionsSection = tractions.section();
  CPPUNIT_ASSERT(!tractionsSection.isNull());

  const PylithScalar t = 0;
  fault.updateStateVars(t, &fields);  
  fault._calcTractionsChange(&tractions, fields.get("disp(t)"));

  int iVertex = 0;
  const PylithScalar tolerance = 1.0e-06;
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
    const PylithScalar* tractionsVertex = tractionsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != tractionsVertex);

    fiberDim = dispSection->getFiberDimension(meshVertex);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    const PylithScalar* dispVertex = dispSection->restrictPoint(meshVertex);
    CPPUNIT_ASSERT(0 != dispVertex);

    const PylithScalar scale = 1.0 / _data->area[iVertex];
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const PylithScalar tractionE = dispVertex[iDim] * scale;
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
// Test splitField().
void
pylith::faults::TestFaultCohesiveKin::testSplitField(void)
{ // testSplitField
  assert(0 != _data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const topology::Field<topology::Mesh>& disp = fields.get("disp(t)");
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);
  const int spaceDim = cs->spaceDim();

  topology::Field<topology::Mesh> splitField(mesh);
  splitField.newSection(disp, spaceDim);
  splitField.splitDefault();
  fault.splitField(&splitField);
  splitField.allocate();
  splitField.zero();

  const ALE::Obj<RealSection>& section = splitField.section();
  CPPUNIT_ASSERT(!section.isNull());
  CPPUNIT_ASSERT_EQUAL(spaceDim+1, section->getNumSpaces());

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  const int numFaultVertices = _data->numFaultVertices;
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    const int fiberDim = section->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(spaceDim, fiberDim);
    bool isLagrangeVertex = false;
    for (int i=0; i < numFaultVertices; ++i)
      if (*v_iter == _data->verticesLagrange[i]) {
	isLagrangeVertex = true;
	break;
      } // if
    if (isLagrangeVertex) {
      for (int fibration=0; fibration < spaceDim; ++fibration)
        CPPUNIT_ASSERT_EQUAL(0, section->getFiberDimension(*v_iter, fibration));
      const int fibrationF = spaceDim;
      CPPUNIT_ASSERT_EQUAL(spaceDim, 
			   section->getFiberDimension(*v_iter, fibrationF));
    } else {
      for (int fibration=0; fibration < spaceDim; ++fibration)
        CPPUNIT_ASSERT_EQUAL(1, section->getFiberDimension(*v_iter, fibration));
      const int fibrationF = spaceDim;
      CPPUNIT_ASSERT_EQUAL(0, section->getFiberDimension(*v_iter, fibrationF));
    } // if/else
  } // for
} // testSplitField

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
  
  //mesh->debug(true); // DEBUGGING
  
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
  
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(_data->label)->size();
  int firstFaultCell      = mesh->sieveMesh()->getIntSection(_data->label)->size();
  if (fault->useLagrangeConstraints()) {
    firstFaultCell += mesh->sieveMesh()->getIntSection(_data->label)->size();
  }
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->eqsrcs(const_cast<const char**>(names), nsrcs, sources, nsrcs);
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);
  
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };
  
  fault->initialize(*mesh, upDir); 
  
  delete[] sources; sources = 0;
  for (int i=0; i < nsrcs; ++i)
    delete[] names[i];
  delete[] names; names = 0;
  
  // Setup fields
  fields->add("residual", "residual");
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("dispIncr adjust", "displacement_adjust");
  fields->solutionName("dispIncr(t->t+dt)");
  
  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  fields->copyLayout("residual");
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveKin::_isLagrangeVertex(const int vertex) const
{ // _isLagrangeVertex
  assert(0 != _data);

  const int numFaultVertices = _data->numFaultVertices;
  bool isFound = false;
  for (int i=0; i < _data->numFaultVertices; ++i)
    if (_data->verticesLagrange[i] == vertex) {
      isFound = true;
      break;
    } // if
  return isFound;
} // _isLagrangeVertex


// End of file 
