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
// Copyright (c) 2010-2013 University of California, Davis
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

  _flipFault = false;

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

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  // Check fault mesh sizes
  CPPUNIT_ASSERT_EQUAL(_data->cellDim, fault.dimension());
  CPPUNIT_ASSERT_EQUAL(_data->numBasis, fault.numCorners());
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, fault.numVertices());
  CPPUNIT_ASSERT_EQUAL(_data->numCohesiveCells, fault.numCells());

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::SubMeshIS subpointIS(*fault._faultMesh);
  const PetscInt numPoints = subpointIS.size();
  const PetscInt* points = subpointIS.points();CPPUNIT_ASSERT(points);

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;

    PetscErrorCode err = PetscFindInt(_data->verticesNegative[v-vStart], numPoints, points, &faultPoint);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, _data->verticesFault[v-vStart]);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, vEnd-vStart);

  // Check cohesive vertex info
  const int numFaultVertices = _data->numFaultVertices;
  CPPUNIT_ASSERT_EQUAL(numFaultVertices, int(fault._cohesiveVertices.size()));
  for (int i=0; i < numFaultVertices; ++i) {
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[i], fault._cohesiveVertices[i].fault);
    CPPUNIT_ASSERT_EQUAL(_data->verticesLagrange[i], fault._cohesiveVertices[i].lagrange);
    CPPUNIT_ASSERT_EQUAL(_data->verticesNegative[i], fault._cohesiveVertices[i].negative);
    CPPUNIT_ASSERT_EQUAL(_data->verticesPositive[i], fault._cohesiveVertices[i].positive);
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

  const int spaceDim = _data->spaceDim;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  CPPUNIT_ASSERT(_data->fieldT);
  topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
  PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = dispVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d] / _data->lengthScale;
    } // for
  } // for
  
  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);
  topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);

  //residual.view("RESIDUAL"); // DEBUGGING

  // Check values
  CPPUNIT_ASSERT(_data->residual);
  const PylithScalar* valsE = _data->residual;
  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);
      
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  const PylithScalar residualScale = _data->lengthScale * pow(_data->lengthScale, spaceDim-1);
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = residualVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, residualVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar valE = valsE[iVertex*spaceDim+d];
      if (fabs(valE) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valE*residualScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+d]*residualScale, tolerance);
    } // for
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

  const int spaceDim = _data->spaceDim;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  CPPUNIT_ASSERT(_data->fieldT);
  topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
  PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = dispVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
    } // for
  } // for
  
  const PylithScalar t = 2.134;
  topology::Jacobian jacobian(fields.solution());
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
  jacobian.assemble("final_assembly");

  //MatView(jacobian.matrix(), PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  CPPUNIT_ASSERT(_data->jacobian);
  const PylithScalar* valsE = _data->jacobian;
  const PetscInt nrowsE = verticesStratum.size() * _data->spaceDim;
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
      if (fabs(valE-vals[index]) > tolerance)
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
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const PylithScalar t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
  jacobian.complete();

  //jacobian.view("JACOBIAN"); // DEBUGGING

  // Only check Lagrange multiplier values
  const int spaceDim = _data->spaceDim;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();

  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  PetscErrorCode  err;
  err = DMPlexGetHybridBounds(dmMesh, &cStart, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();CPPUNIT_ASSERT(jacobianArray);
  
  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone = NULL;
    PetscInt coneSize, p;

    err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    // Check Lagrange multiplier dofs
    //   For depth = 1, we have a prism and use the last third
    coneSize /= 3;
    //   For depth > 1, we take the edges
    for (p = 2*coneSize; p < 3*coneSize; ++p) {
      const PetscInt v = cone[p];

      const PetscInt off = jacobianVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, jacobianVisitor.sectionDof(v));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        const PylithScalar valE = _data->jacobianLumped[(v-vStart)*spaceDim+d];
#if 0 // debugging
        std::cout << "vertex: " << *v_iter << ", iDim: " << iDim
                  << ", valE: " << valE
                  << ", val: " << vals[iDim]
                  << std::endl;
#endif
        if (fabs(valE) > 1.0)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobianArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, jacobianArray[off+d], tolerance);
      } // for
    } // for
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
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  { // setup disp
    CPPUNIT_ASSERT(_data->fieldT);
    topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
    PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = dispVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(v));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d] / _data->lengthScale;
      } // for
    } // for
  } // setup disp

  // compute residual so that slip and residual are setup
  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);
  topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);
  residual.complete();

  { // setup disp increment
    CPPUNIT_ASSERT(_data->fieldIncr);
    topology::VecVisitorMesh dispIncrVisitor(fields.get("dispIncr(t->t+dt)"));
    PetscScalar* dispIncrArray = dispIncrVisitor.localArray();CPPUNIT_ASSERT(dispIncrArray);
    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = dispIncrVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispIncrVisitor.sectionDof(v));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispIncrArray[off+d] = _data->fieldIncr[iVertex*spaceDim+d] / _data->lengthScale;
      } // for
    } // for
  } // setup disp increment

  // Set Jacobian values
  topology::Field jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();
  { // setup jacobian
    const PylithScalar jacobianScale = pow(_data->lengthScale, spaceDim-1);
    topology::VecVisitorMesh jacobianVisitor(jacobian);
    PetscScalar* jacobianArray = jacobianVisitor.localArray();CPPUNIT_ASSERT(jacobianArray);
    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = jacobianVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, jacobianVisitor.sectionDof(v));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        jacobianArray[off+d] = _data->jacobianLumped[iVertex*spaceDim+d] / jacobianScale;
      } // for
    } // for
  } // setup jacobian
  jacobian.complete();

  topology::Field& solution = fields.get("dispIncr(t->t+dt)");
  fault.adjustSolnLumped(&fields, t, jacobian);
  const topology::Field& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  //solution.view("SOLUTION AFTER ADJUSTMENT"); // DEBUGGING

  topology::VecVisitorMesh solutionVisitor(solution);
  const PetscScalar* solutionArray = solutionVisitor.localArray();CPPUNIT_ASSERT(solutionArray);

  CPPUNIT_ASSERT(_data->fieldIncrAdjusted);
  const PylithScalar* solutionE = _data->fieldIncrAdjusted;
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  for(PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = solutionVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, solutionVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d, ++index) {
      if (0.0 != solutionE[index])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, solutionArray[off+d]/solutionE[index]*_data->lengthScale, tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(solutionE[index], solutionArray[off+d]*_data->lengthScale, tolerance);
    } // for
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

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const int spaceDim = _data->spaceDim;
  { // setup disp
    topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
    PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = dispVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(v));
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
      } // for
    } // for
  } // setup disp

  CPPUNIT_ASSERT(fault._faultMesh);
  topology::Field tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();

  const PylithScalar t = 0;
  fault.updateStateVars(t, &fields);  
  fault._calcTractionsChange(&tractions, fields.get("disp(t)"));

  PetscDM faultDMMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);

  topology::SubMeshIS subpointIS(*fault._faultMesh);
  const PetscInt numPoints = subpointIS.size();
  const PetscInt* points = subpointIS.points();CPPUNIT_ASSERT(points);

  topology::VecVisitorMesh tractionVisitor(tractions);
  const PetscScalar* tractionArray = tractionVisitor.localArray();CPPUNIT_ASSERT(tractionArray);

  topology::VecVisitorMesh dispVisitor(fields.get("disp(t)"));
  const PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);

  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  PetscErrorCode err = DMPlexGetHybridBounds(dmMesh, &cStart, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);

  topology::Stratum fverticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt fvStart = fverticesStratum.begin();

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone = NULL;
    PetscInt coneSize = 0, p = 0;

    err = DMPlexGetConeSize(dmMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetCone(dmMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    // Check Lagrange multiplier dofs
    //   For depth = 1, we have a prism and use the last third
    CPPUNIT_ASSERT_EQUAL(0, coneSize % 3);
    coneSize /= 3;
    //   For depth > 1, we take the edges
    for (p = 2*coneSize; p < 3*coneSize; ++p) {
      const PetscInt v_lagrange = cone[p];
      const PetscInt v_negative = cone[p-2*coneSize];
      PetscInt v_fault;

      err = PetscFindInt(v_negative, numPoints, points, &v_fault);PYLITH_CHECK_ERROR(err);CPPUNIT_ASSERT(v_fault >= 0);
      const PetscInt toff = tractionVisitor.sectionOffset(v_fault);
      CPPUNIT_ASSERT_EQUAL(spaceDim, tractionVisitor.sectionDof(v_fault));

      const PetscInt doff = dispVisitor.sectionOffset(v_lagrange);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(v_lagrange));
      
      const PylithScalar* orientationVertex = &_data->orientation[(v_fault-fvStart)*spaceDim*spaceDim];

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
  } // for

  PYLITH_METHOD_END;
} // testCalcTractionsChange

// ----------------------------------------------------------------------
// Test splitField().
void
pylith::faults::TestFaultCohesiveKin::testSplitField(void)
{ // testSplitField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const topology::Field& disp = fields.get("disp(t)");
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field splitField(mesh);
  splitField.addField("displacement", spaceDim);
  splitField.addField("multipliers", spaceDim);
  splitField.setupFields();
  splitField.newSection(disp, spaceDim);
  fault.splitField(&splitField);
  splitField.allocate();
  splitField.zero();

  PetscSection section = splitField.petscSection();CPPUNIT_ASSERT(section);
  PetscVec vec = splitField.localVector();CPPUNIT_ASSERT(vec);
  PetscScalar *array = NULL;
  PetscInt numFields, numComp;
  PetscErrorCode err;
  err = PetscSectionGetNumFields(section, &numFields);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(2, numFields);
  err = PetscSectionGetFieldComponents(section, 0, &numComp);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(spaceDim, numComp);
  err = PetscSectionGetFieldComponents(section, 1, &numComp);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(spaceDim, numComp);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const int numFaultVertices = _data->numFaultVertices;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, fdof;

    err = PetscSectionGetDof(section, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(section, v, &off);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    bool isLagrangeVertex = false;
    for(int i = 0; i < numFaultVertices; ++i)
      if (v == _data->verticesLagrange[i]) {
        isLagrangeVertex = true;
        break;
      } // if
    if (isLagrangeVertex) {
      err = PetscSectionGetFieldDof(section, v, 0, &fdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(0, fdof);
      err = PetscSectionGetFieldDof(section, v, 1, &fdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fdof);
    } else {
      err = PetscSectionGetFieldDof(section, v, 0, &fdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fdof);
      err = PetscSectionGetFieldDof(section, v, 1, &fdof);PYLITH_CHECK_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(0, fdof);
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testSplitField

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
  fault->quadrature(_quadrature);
  fault->eqsrcs(const_cast<const char**>(names), nsrcs, sources, nsrcs);

  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = 0;
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err = DMPlexGetStratumSize(dmMesh, _data->label, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
  PetscInt firstFaultCell = firstLagrangeVertex + firstLagrangeVertex;
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);
  
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
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  fields->copyLayout("residual");

  fault->verifyConfiguration(*mesh);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveKin::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  const int numFaultVertices = _data->numFaultVertices;
  bool isFound = false;
  for (int i=0; i < _data->numFaultVertices; ++i)
    if (_data->verticesLagrange[i] == vertex) {
      isFound = true;
      break;
    } // if
  PYLITH_METHOD_RETURN(isFound);
} // _isConstraintVertex


// End of file 
