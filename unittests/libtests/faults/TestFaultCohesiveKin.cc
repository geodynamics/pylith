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
// Copyright (c) 2010-2015 University of California, Davis
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
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES SubMeshIS
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKin );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKin::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _quadrature = new feassemble::Quadrature();CPPUNIT_ASSERT(_quadrature);
  const int nsrcs = 1;
  _eqsrcs.resize(nsrcs);
  _eqsrcs[0] = new EqKinSrc();
  _slipfns.resize(nsrcs);
  _slipfns[0] = new BruneSlipFn();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveKin::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  int nsrcs = _eqsrcs.size();
  for (int i=0; i < nsrcs; ++i)
    delete _eqsrcs[i];
  nsrcs = _slipfns.size();
  for (int i=0; i < nsrcs; ++i)
    delete _slipfns[i];

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveKin::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  FaultCohesiveKin fault;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test eqsrc().
void
pylith::faults::TestFaultCohesiveKin::testEqsrc(void)
{ // testEqsrc
  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_END;
} // testEqsrc

// ----------------------------------------------------------------------
/// Test needNewJacobian()
void
pylith::faults::TestFaultCohesiveKin::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  PYLITH_METHOD_BEGIN;

  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault.needNewJacobian());
  fault._needNewJacobian = false;
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  PYLITH_METHOD_END;
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test useLagrangeConstraints().
void
pylith::faults::TestFaultCohesiveKin::testUseLagrangeConstraints(void)
{ // testUseLagrangeConstraints
  PYLITH_METHOD_BEGIN;

  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault.useLagrangeConstraints());

  PYLITH_METHOD_END;
} // testUseLagrangeConstraints

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveKin::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  PetscErrorCode err;

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

#if 0 // DEBUGGING
  mesh.view("::ascii_info_detail");
  fault._faultMesh->view("::ascii_info_detail");
#endif

  // Check fault mesh sizes
  CPPUNIT_ASSERT_EQUAL(_data->cellDim, fault.dimension());
  CPPUNIT_ASSERT_EQUAL(_data->numBasis, fault.numCorners());
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, fault.numVertices());
  CPPUNIT_ASSERT_EQUAL(_data->numCohesiveCells, fault.numCells());

  topology::SubMeshIS subpointIS(*fault._faultMesh);
  const PetscInt numPoints = subpointIS.size();
  const PetscInt* points = subpointIS.points();CPPUNIT_ASSERT(points);

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;

    err = PetscFindInt(_data->verticesNegative[v-vStart], numPoints, points, &faultPoint);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, _data->verticesFault[v-vStart]);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, vEnd-vStart);

  // Check cohesive vertex info; permit different order of vertices.
  PetscDMLabel clamped = NULL;
  err = DMGetLabel(dmMesh, "clamped", &clamped);PYLITH_CHECK_ERROR(err);

  const int numFaultVertices = _data->numFaultVertices;
  CPPUNIT_ASSERT_EQUAL(numFaultVertices, int(fault._cohesiveVertices.size()));
  int_array verticesFaultSorted(_data->verticesFault, numFaultVertices);
  int* sortedBegin = &verticesFaultSorted[0];
  int* sortedEnd = &verticesFaultSorted[numFaultVertices];
  std::sort(sortedBegin, sortedEnd);
  for (int i=0; i < numFaultVertices; ++i) {
    const PetscInt sign = (fault._cohesiveVertices[i].fault < 0) ? -1 : 1;
    const PetscInt v_fault = sign*fault._cohesiveVertices[i].fault;
    const int* iter = std::lower_bound(sortedBegin, sortedEnd, v_fault);
    CPPUNIT_ASSERT(iter != sortedEnd);
    const int index = iter - sortedBegin;

    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[index], sign*fault._cohesiveVertices[i].fault);
    CPPUNIT_ASSERT_EQUAL(_data->edgesLagrange[index], sign*fault._cohesiveVertices[i].lagrange);
    CPPUNIT_ASSERT_EQUAL(_data->verticesNegative[index], sign*fault._cohesiveVertices[i].negative);
    CPPUNIT_ASSERT_EQUAL(_data->verticesPositive[index], sign*fault._cohesiveVertices[i].positive);
  } // for

  // Check cohesive cell info
  const int numCohesiveCells = _data->numCohesiveCells;
  CPPUNIT_ASSERT_EQUAL(numCohesiveCells, int(fault._cohesiveToFault.size()));
  std::map<PetscInt,PetscInt>::iterator m_iterator = fault._cohesiveToFault.begin();
  std::map<PetscInt,PetscInt>::const_iterator mapEnd = fault._cohesiveToFault.end();
  for (int i=0; i < numCohesiveCells; ++i, ++m_iterator) {
    CPPUNIT_ASSERT(mapEnd != m_iterator);
    CPPUNIT_ASSERT_EQUAL(_data->cellMappingFault[i], m_iterator->second);
    CPPUNIT_ASSERT_EQUAL(_data->cellMappingCohesive[i], m_iterator->first);
  } // for

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING
  topology::VecVisitorMesh orientationVisitor(fault._fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();CPPUNIT_ASSERT(orientationArray);

  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = orientationVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(orientationSize, orientationVisitor.sectionDof(v));

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < orientationSize; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[iVertex*orientationSize+d], orientationArray[off+d], tolerance);
    } // for
  } // for

  // Check area
  topology::VecVisitorMesh areaVisitor(fault._fields->get("area"));
  const PetscScalar* areaArray = areaVisitor.localArray();CPPUNIT_ASSERT(areaArray);
  const PylithScalar areaScale = pow(_data->lengthScale, spaceDim-1);
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = areaVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, areaVisitor.sectionDof(v));
    const PylithScalar tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->area[iVertex], areaArray[off]*areaScale, tolerance);
  } // for

  PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  CPPUNIT_ASSERT(_data->fieldT);
  _fieldSetValues(&fields.get("disp(t)"), _data->fieldT, _data->lengthScale);
  
  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);
  topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);

  //residual.view("RESIDUAL"); // DEBUGGING

  // Check values
  CPPUNIT_ASSERT(_data->residual);
  const PylithScalar* valsE = _data->residual;

  PetscInt pStart, pEnd;
  PetscErrorCode err = PetscSectionGetChart(residual.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);
      
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  const int spaceDim = _data->spaceDim;
  const PylithScalar residualScale = _data->lengthScale * pow(_data->lengthScale, spaceDim-1);
  for (PetscInt p = pStart, iPoint = 0; p < pEnd; ++p) {
    if (residualVisitor.sectionDof(p) > 0) {
      const PetscInt off = residualVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, residualVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	const PylithScalar valE = valsE[iPoint*spaceDim+d];
	if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valE*residualScale, tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+d]*residualScale, tolerance);
      } // for
      ++iPoint;
    } // if
  } // for


  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  CPPUNIT_ASSERT(_data->fieldT);
  _fieldSetValues(&fields.get("disp(t)"), _data->fieldT);
  
  const PylithScalar t = 2.134;
  topology::Jacobian jacobian(fields.solution());
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
  jacobian.assemble("final_assembly");

  //MatView(jacobian.matrix(), PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);

  PetscInt numClampedVertices = 0;
  PetscErrorCode err = DMGetLabelSize(fault._faultMesh->dmMesh(), "clamped", &numClampedVertices);PYLITH_CHECK_ERROR(err);

  CPPUNIT_ASSERT(_data->jacobian);
  const int numDOF = verticesStratum.size() + _data->numFaultVertices - numClampedVertices;
  const int spaceDim = _data->spaceDim;
  const PylithScalar* valsE = _data->jacobian;
  const PetscInt nrowsE = numDOF * _data->spaceDim;
  const PetscInt ncolsE = nrowsE;

  PetscMat jacobianMat = jacobian.matrix();
  int nrows = 0;
  int ncols = 0;
  MatGetSize(jacobianMat, &nrows, &ncols);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  MatConvert(jacobianMat, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);

  scalar_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar jacobianScale = pow(_data->lengthScale, spaceDim-1);
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      const PylithScalar valE = valsE[index];
#if 0 // DEBUGGING
      if (fabs(valE-vals[index]*jacobianScale) > tolerance)
	std::cout << "ERROR: iRow: " << iRow << ", iCol: " << iCol
		  << "valE: " << valE
		  << ", val: " << vals[index]
		  << std::endl;
#endif // DEBUGGING
      if (fabs(valE) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valE*jacobianScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[index]*jacobianScale, tolerance);
    } // for
  MatDestroy(&jDense);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  PYLITH_METHOD_END;
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobian() with lumped Jacobian.
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->jacobianLumped);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  topology::Field jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.cloneSection(fields.get("residual"));

  const PylithScalar t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
  jacobian.complete();

  //jacobian.view("JACOBIAN"); // DEBUGGING

  const int spaceDim = _data->spaceDim;

  PetscInt pStart, pEnd;
  PetscErrorCode err = PetscSectionGetChart(jacobian.localSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();CPPUNIT_ASSERT(jacobianArray);
  
  const PylithScalar tolerance = 1.0e-06;
  for (PetscInt p = pStart, iPoint = 0; p < pEnd; ++p) {
    if (jacobianVisitor.sectionDof(p) > 0) {
      const PetscInt off = jacobianVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, jacobianVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        const PylithScalar valE = (p >= vEnd) ? 1.0 : 0.0; // Lagrange multiplier DOF should be 1.0, else 0.0.
	if (fabs(valE) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobianArray[off+d]/valE, tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, jacobianArray[off+d], tolerance);
      } // for
      ++iPoint;
    } // if
  } // for

  PYLITH_METHOD_END;
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
// Test adjustSolnLumped().
void
pylith::faults::TestFaultCohesiveKin::testAdjustSolnLumped(void)
{ // testAdjustSolnLumped
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->jacobianLumped);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;

  _fieldSetValues(&fields.get("disp(t)"), _data->fieldT, _data->lengthScale);

  // compute residual so that slip and residual are setup
  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);
  topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);
  residual.complete();

  // Setup disp increment
  _fieldSetValues(&fields.get("dispIncr(t->t+dt)"), _data->fieldIncr, _data->lengthScale);

  // Set Jacobian values
  const PylithScalar jacobianScale = pow(_data->lengthScale, spaceDim-1);
  topology::Field jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.cloneSection(fields.get("residual"));
  _fieldSetValues(&jacobian, _data->jacobianLumped, jacobianScale);
  jacobian.complete();

  topology::Field& solution = fields.get("dispIncr(t->t+dt)");
  fault.adjustSolnLumped(&fields, t, jacobian);
  const topology::Field& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  //solution.view("SOLUTION AFTER ADJUSTMENT"); // DEBUGGING

  PetscInt pStart, pEnd;
  PetscErrorCode err = PetscSectionGetChart(jacobian.localSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
  topology::VecVisitorMesh solutionVisitor(solution);
  const PetscScalar* solutionArray = solutionVisitor.localArray();CPPUNIT_ASSERT(solutionArray);

  CPPUNIT_ASSERT(_data->fieldIncrAdjusted);
  const PylithScalar* solutionE = _data->fieldIncrAdjusted;
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;

  for (PetscInt p = pStart, iPoint = 0; p < pEnd; ++p) {
    if (solutionVisitor.sectionDof(p) > 0) {
      const PetscInt off = solutionVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, solutionVisitor.sectionDof(p));
      for (PetscInt d = 0; d < spaceDim; ++d) {
	if (0.0 != solutionE[iPoint*spaceDim+d])
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, solutionArray[off+d]/solutionE[iPoint*spaceDim+d]*_data->lengthScale, tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(solutionE[iPoint*spaceDim+d], solutionArray[off+d]*_data->lengthScale, tolerance);
      } // for
      ++iPoint;
    } // if
  } // for

  PYLITH_METHOD_END;
} // testAdjustSolnLumped

// ----------------------------------------------------------------------
// Test calcTractionsChange().
void
pylith::faults::TestFaultCohesiveKin::testCalcTractionsChange(void)
{ // testCalcTractionsChange
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  // Setup dispT
  _fieldSetValues(&fields.get("disp(t)"), _data->fieldT);

  const int spaceDim = _data->spaceDim;

  CPPUNIT_ASSERT(fault._faultMesh);
  topology::Field tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();

  const PylithScalar t = 0;
  fault.updateStateVars(t, &fields);  
  fault._calcTractionsChange(&tractions, fields.get("disp(t)"));

  topology::VecVisitorMesh tractionVisitor(tractions);
  const PetscScalar* tractionArray = tractionVisitor.localArray();CPPUNIT_ASSERT(tractionArray);

  topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
  const PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);

  const int numFaultVertices = _data->numFaultVertices;
  CPPUNIT_ASSERT_EQUAL(numFaultVertices, int(fault._cohesiveVertices.size()));
  int_array verticesFaultSorted(_data->verticesFault, numFaultVertices);
  int* sortedBegin = &verticesFaultSorted[0];
  int* sortedEnd = &verticesFaultSorted[numFaultVertices];
  std::sort(sortedBegin, sortedEnd);
  for (int i=0; i < numFaultVertices; ++i) {
    const int* iter = std::lower_bound(sortedBegin, sortedEnd, fault._cohesiveVertices[i].fault);
    CPPUNIT_ASSERT(iter != sortedEnd);
    const int index = iter - sortedBegin;

    const PetscInt v_fault = fault._cohesiveVertices[i].fault;
    const PetscInt e_lagrange = fault._cohesiveVertices[i].lagrange;
    if (e_lagrange < 0) { // skip clamped edges
      continue;
    } // if

    const PetscInt toff = tractionVisitor.sectionOffset(v_fault);
    CPPUNIT_ASSERT_EQUAL(spaceDim, tractionVisitor.sectionDof(v_fault));

    const PetscInt doff = dispVisitor.sectionOffset(e_lagrange);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(e_lagrange));
      
    const PylithScalar* orientationVertex = &_data->orientation[index*spaceDim*spaceDim];

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < spaceDim; ++d) {
      PylithScalar tractionE = 0.0;
      for(PetscInt e = 0; e < spaceDim; ++e)
	tractionE += orientationVertex[d*spaceDim+e] * dispArray[doff+e];
      if (tractionE > 1.0) 
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tractionArray[toff+d]/tractionE, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, tractionArray[toff+d], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testCalcTractionsChange


// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesiveKin::_fieldSetValues(topology::Field* field,
						      const PylithScalar* data,
						      const PylithScalar scale)
{ // _fieldSetValues
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(field);
  CPPUNIT_ASSERT(data);
  CPPUNIT_ASSERT(_data);

  const int spaceDim = _data->spaceDim;

  PetscErrorCode err;
  PetscInt pStart, pEnd;
  err = PetscSectionGetChart(field->localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);

  topology::VecVisitorMesh fieldVisitor(*field);
  PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
  for (PetscInt p = pStart, iPoint = 0; p < pEnd; ++p) {
    if (fieldVisitor.sectionDof(p) > 0) {
      const PetscInt off = fieldVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fieldVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	fieldArray[off+d] = data[iPoint*spaceDim+d] / scale;
      } // for
      ++iPoint;
    } // if
  } // for
} // _fieldSetValues


// ----------------------------------------------------------------------
// Initialize FaultCohesiveKin interface condition.
void
pylith::faults::TestFaultCohesiveKin::_initialize(topology::Mesh* const mesh,
						  FaultCohesiveKin* const fault,
						  topology::SolutionFields* const fields) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  
  // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Set scales
  // Most test data is insensitive to the scales because we set the fields directly.
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);
  
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

  fault->id(_data->id);
  fault->label(_data->label);
  if (_data->edge) {
    fault->edge(_data->edge);
  } // if
  fault->quadrature(_quadrature);
  fault->eqsrcs(const_cast<const char**>(names), nsrcs, sources, nsrcs);

  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = 0;
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err = DMGetStratumSize(dmMesh, _data->label, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
  PetscInt firstFaultCell = firstLagrangeVertex + firstLagrangeVertex;
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  fault->normalizer(normalizer);
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
  topology::Field& residual = fields->get("residual");
  residual.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR);
  residual.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);
  residual.subfieldsSetup();
  residual.setupSolnChart();
  residual.setupSolnDof(spaceDim);
  fault->setupSolnDof(&residual);
  residual.allocate();
  residual.zeroAll();

  fields->copyLayout("residual");

  fault->verifyConfiguration(*mesh);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Determine if point is a Lagrange multiplier constraint point.
bool
pylith::faults::TestFaultCohesiveKin::_isConstraintEdge(const int point) const
{ // _isConstraintEdge
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  const int numFaultVertices = _data->numFaultVertices;
  bool isFound = false;
  for (int i=0; i < _data->numFaultVertices; ++i)
    if (_data->edgesLagrange[i] == point) {
      isFound = true;
      break;
    } // if
  PYLITH_METHOD_RETURN(isFound);
} // _isConstraintEdge


// End of file 
