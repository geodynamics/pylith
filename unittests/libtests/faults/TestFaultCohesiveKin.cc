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
// Copyright (c) 2010-2012 University of California, Davis
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
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  // Check fault mesh sizes
  CPPUNIT_ASSERT_EQUAL(_data->cellDim, fault.dimension());
  CPPUNIT_ASSERT_EQUAL(_data->numBasis, fault.coneSize());
  CPPUNIT_ASSERT_EQUAL(_data->numFaultVertices, fault.numVertices());
  CPPUNIT_ASSERT_EQUAL(_data->numCohesiveCells, fault.numCells());

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscIS subpointIS;
  const PetscInt *points;
  PetscInt vStart, vEnd, numPoints;
  PetscErrorCode  err;

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(dmMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(subpointIS);
  err = ISGetSize(subpointIS, &numPoints);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;

    err = PetscFindInt(_data->verticesLagrange[v-vStart], numPoints, points, &faultPoint);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, _data->verticesFault[v-vStart]);
  } // for
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
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
  PetscSection orientationSection = fault._fields->get("orientation").petscSection();
  Vec          orientationVec     = fault._fields->get("orientation").localVector();
  PetscScalar *orientationArray;
  CPPUNIT_ASSERT(orientationSection);CPPUNIT_ASSERT(orientationVec);
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  int iVertex = 0;
  err = VecGetArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt dof, off;

    err = PetscSectionGetDof(orientationSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(orientationSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(orientationSize, dof);

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < orientationSize; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[iVertex*orientationSize+d], orientationArray[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(orientationVec, &orientationArray);CHECK_PETSC_ERROR(err);

  // Check area
  PetscSection areaSection = fault._fields->get("area").petscSection();
  Vec          areaVec     = fault._fields->get("area").localVector();
  PetscScalar *areaArray;
  CPPUNIT_ASSERT(areaSection);CPPUNIT_ASSERT(areaVec);
  iVertex = 0;
  err = VecGetArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt dof, off;

    err = PetscSectionGetDof(areaSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(areaSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(1, dof);

    const PylithScalar tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->area[iVertex], areaArray[off], tolerance);
  } // for
  err = VecRestoreArray(areaVec, &areaArray);CHECK_PETSC_ERROR(err);
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);
  CPPUNIT_ASSERT(_data->residual);
  CPPUNIT_ASSERT(_data->residualIncr);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  PetscScalar *residualArray;
  CPPUNIT_ASSERT(residualSection);CPPUNIT_ASSERT(residualVec);

  PetscSection dispSection = fields.get("disp(t)").petscSection();
  Vec          dispVec     = fields.get("disp(t)").localVector();
  PetscScalar *dispArray;
  CPPUNIT_ASSERT(dispSection);CPPUNIT_ASSERT(dispVec);

  DM              dmMesh = mesh.dmMesh();
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  int iVertex = 0;
  err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt dof, off;

    err = PetscSectionGetDof(dispSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    for(PetscInt d = 0; d < dof; ++d) {
      dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
    }
  }
  err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  
  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);

  { // Integrate residual with disp (as opposed to disp increment).
    fault.integrateResidual(residual, t, &fields);

    //residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const PylithScalar* valsE = _data->residual;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(residualSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(residualSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);
      
      for(PetscInt d = 0; d < fiberDimE; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  } // Integrate residual with disp (as opposed to disp increment).

  residual.zero();
  { // Integrate residual with disp increment.
    fault.integrateResidual(residual, t, &fields);

    //residual.view("RESIDUAL"); // DEBUGGING

    // Check values
    const PylithScalar* valsE = _data->residualIncr;
    iVertex = 0;
    const int fiberDimE = spaceDim;
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(residualSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(residualSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);
      
      for(PetscInt d = 0; d < fiberDimE; ++d) {
        const PylithScalar valE = valsE[iVertex*spaceDim+d];
        if (fabs(valE) > tolerance)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valE, tolerance);
        else
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+d], tolerance);
      } // for
    } // for
    err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  } // Integrate residual with disp increment.
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);
  CPPUNIT_ASSERT(_data->jacobian);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  PetscSection dispSection = fields.get("disp(t)").petscSection();
  Vec          dispVec     = fields.get("disp(t)").localVector();
  PetscScalar *dispArray;
  CPPUNIT_ASSERT(dispSection);CPPUNIT_ASSERT(dispVec);

  DM              dmMesh = mesh.dmMesh();
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  int iVertex = 0;
  err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt dof, off;

    err = PetscSectionGetDof(dispSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
    }
  }
  err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);

  topology::Jacobian jacobian(fields.solution());

  const PylithScalar t = 2.134;
  fault.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  jacobian.assemble("final_assembly");

  //MatView(jacobian.matrix(), PETSC_VIEWER_STDOUT_WORLD); // DEBUGGING

  const PylithScalar* valsE = _data->jacobian;
  PetscInt nrowsE, ncolsE;
  err = PetscSectionGetStorageSize(dispSection, &nrowsE);CHECK_PETSC_ERROR(err);
  ncolsE = nrowsE;

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
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobian() with lumped Jacobian.
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->jacobianLumped);

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

  //jacobian.view("JACOBIAN"); // DEBUGGING

  PetscSection jacobianSection = jacobian.petscSection();
  Vec          jacobianVec     = jacobian.localVector();
  PetscScalar *jacobianArray;
  CPPUNIT_ASSERT(jacobianSection);CPPUNIT_ASSERT(jacobianVec);

  // Only check Lagrange multiplier values

  const int spaceDim = _data->spaceDim;

  DM              dmMesh = mesh.dmMesh();
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  DM              faultDMMesh = fault._faultMesh->dmMesh();
  IS              subpointIS;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(subpointIS);

  const PylithScalar tolerance = 1.0e-06;
  int iVertex = 0;
  const PetscInt *points;
  PetscInt        numPoints;

  err = ISGetSize(subpointIS, &numPoints);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    PetscInt faultPoint;

    err = PetscFindInt(v, numPoints, points, &faultPoint);CHECK_PETSC_ERROR(err);
    if (faultPoint < 0) // only check Lagrange multiplier values
      continue;
    PetscInt dof, off;

    err = PetscSectionGetDof(jacobianSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(jacobianSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    
    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar valE = _data->jacobianLumped[iVertex*spaceDim+d];
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
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
// Test adjustSolnLumped().
void
pylith::faults::TestFaultCohesiveKin::testAdjustSolnLumped(void)
{ // testAdjustSolnLumped
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);
  CPPUNIT_ASSERT(_data->fieldIncr);
  CPPUNIT_ASSERT(_data->fieldIncrAdjusted);
  CPPUNIT_ASSERT(_data->jacobianLumped);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const int spaceDim = _data->spaceDim;
  DM              dmMesh = mesh.dmMesh();
  PetscInt        vStart, vEnd;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  { // setup disp
    PetscSection dispSection = fields.get("disp(t)").petscSection();
    Vec          dispVec     = fields.get("disp(t)").localVector();
    PetscScalar *dispArray;
    CPPUNIT_ASSERT(dispSection);CPPUNIT_ASSERT(dispVec);
    int iVertex = 0;
    err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(dispSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
      for(PetscInt d = 0; d < dof; ++d) {
        dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
      }
    } // for
    err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  } // setup disp

  // compute residual so that slip and residual are setup
  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  topology::Field<topology::Mesh>& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);
  residual.complete();

  { // setup disp increment
    PetscSection dispSection = fields.get("dispIncr(t->t+dt)").petscSection();
    Vec          dispVec     = fields.get("dispIncr(t->t+dt)").localVector();
    PetscScalar *dispArray;
    CPPUNIT_ASSERT(dispSection);CPPUNIT_ASSERT(dispVec);
    int iVertex = 0;
    err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(dispSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
      for(PetscInt d = 0; d < dof; ++d) {
        dispArray[off+d] = _data->fieldIncr[iVertex*spaceDim+d];
      }
    } // for
    err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  } // setup disp increment

  // Set Jacobian values
  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();
  { // setup jacobian
    PetscSection jacobianSection = jacobian.petscSection();
    Vec          jacobianVec     = jacobian.localVector();
    PetscScalar *jacobianArray;
    CPPUNIT_ASSERT(jacobianSection);CPPUNIT_ASSERT(jacobianVec);
    int iVertex = 0;
    err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(jacobianSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(jacobianSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
      for(PetscInt d = 0; d < dof; ++d) {
        jacobianArray[off+d] = _data->jacobianLumped[iVertex*spaceDim+d];
      }
    } // for
    err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  } // setup jacobian
  jacobian.complete();

  topology::Field<topology::Mesh>& solution = fields.get("dispIncr(t->t+dt)");
  fault.adjustSolnLumped(&fields, t, jacobian);
  const topology::Field<topology::Mesh>& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  //solution.view("SOLUTION AFTER ADJUSTMENT"); // DEBUGGING

  PetscSection solutionSection = solution.petscSection();
  Vec          solutionVec     = solution.localVector();
  PetscScalar *solutionArray;
  CPPUNIT_ASSERT(solutionSection);CPPUNIT_ASSERT(solutionVec);

  int i = 0;
  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  const PylithScalar* solutionE = _data->fieldIncrAdjusted;
  err = VecGetArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(solutionSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(solutionSection, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);

    for(PetscInt d = 0; d < dof; ++d, ++i) {
      if (0.0 != solutionE[i])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, solutionArray[off+d]/solutionE[i], tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(solutionE[i], solutionArray[off+d], tolerance);
    } // for
  } // for
  err = VecRestoreArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
} // testAdjustSolnLumped

// ----------------------------------------------------------------------
// Test calcTractionsChange().
void
pylith::faults::TestFaultCohesiveKin::testCalcTractionsChange(void)
{ // testCalcTractionsChange
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->fieldT);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  
  const int spaceDim = _data->spaceDim;
  PetscSection dispSection = fields.get("disp(t)").petscSection();
  Vec          dispVec     = fields.get("disp(t)").localVector();
  PetscScalar *dispArray;
  PetscErrorCode err;
  CPPUNIT_ASSERT(dispSection);CPPUNIT_ASSERT(dispVec);
  { // setup disp
    DM              dmMesh = mesh.dmMesh();
    PetscInt        vStart, vEnd;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    int iVertex = 0;
    err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
      PetscInt dof, off;

      err = PetscSectionGetDof(dispSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispSection, v, &off);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
      for(PetscInt d = 0; d < dof; ++d) {
        dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d];
      }
    }
    err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  } // setup disp

  CPPUNIT_ASSERT(0 != fault._faultMesh);
  topology::Field<topology::SubMesh> tractions(*fault._faultMesh);
  tractions.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  tractions.allocate();
  tractions.zero();
  PetscSection tractionsSection = tractions.petscSection();
  Vec          tractionsVec     = tractions.localVector();
  PetscScalar *tractionsArray;
  CPPUNIT_ASSERT(tractionsSection);CPPUNIT_ASSERT(tractionsVec);

  const PylithScalar t = 0;
  fault.updateStateVars(t, &fields);  
  fault._calcTractionsChange(&tractions, fields.get("disp(t)"));

  int iVertex = 0;
  const PylithScalar tolerance = 1.0e-06;
  DM              faultDMMesh = fault._faultMesh->dmMesh();
  IS              subpointIS;
  const PetscInt *points;
  PetscInt        numPoints, vStart, vEnd;

  CPPUNIT_ASSERT(faultDMMesh);
  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT(subpointIS);
  err = ISGetSize(subpointIS, &numPoints);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = VecGetArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v, ++iVertex) {
    const PetscInt meshVertex = points[v];
    PetscInt dof, off, ddof, doff;

    err = PetscSectionGetDof(tractionsSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(tractionsSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(dispSection, meshVertex, &ddof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispSection, meshVertex, &doff);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    CPPUNIT_ASSERT_EQUAL(spaceDim, ddof);
    const PylithScalar *orientationVertex = &_data->orientation[iVertex*spaceDim*spaceDim];

    for(PetscInt d = 0; d < spaceDim; ++d) {
      PylithScalar tractionE = 0.0;
      for(PetscInt e = 0; e < spaceDim; ++e)
        tractionE += orientationVertex[d*spaceDim+e] * dispArray[doff+e];
      if (tractionE > 1.0) 
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tractionsArray[off+d]/tractionE, tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, tractionsArray[off+d], tolerance);
    } // for
  } // for
  err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(tractionsVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispVec, &dispArray);CHECK_PETSC_ERROR(err);
} // testCalcTractionsChange

// ----------------------------------------------------------------------
// Test splitField().
void
pylith::faults::TestFaultCohesiveKin::testSplitField(void)
{ // testSplitField
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  FaultCohesiveKin fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  const topology::Field<topology::Mesh>& disp = fields.get("disp(t)");
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);
  const int spaceDim = cs->spaceDim();

  topology::Field<topology::Mesh> splitField(mesh);
  splitField.addField("displacement", spaceDim);
  splitField.addField("multipliers", spaceDim);
  splitField.setupFields();
  splitField.newSection(disp, spaceDim);
  fault.splitField(&splitField);
  splitField.allocate();
  splitField.zero();

  PetscSection section = splitField.petscSection();
  Vec          vec     = splitField.localVector();
  PetscScalar *array;
  PetscInt     numFields, numComp;
  PetscErrorCode err;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  err = PetscSectionGetNumFields(section, &numFields);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(2, numFields);
  err = PetscSectionGetFieldComponents(section, 0, &numComp);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(spaceDim, numComp);
  err = PetscSectionGetFieldComponents(section, 1, &numComp);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(spaceDim, numComp);

  DM              dmMesh = mesh.dmMesh();
  PetscInt        vStart, vEnd;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const int numFaultVertices = _data->numFaultVertices;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, fdof;

    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dof);
    bool isLagrangeVertex = false;
    for(int i = 0; i < numFaultVertices; ++i)
      if (v == _data->verticesLagrange[i]) {
        isLagrangeVertex = true;
        break;
      } // if
    if (isLagrangeVertex) {
      err = PetscSectionGetFieldDof(section, v, 0, &fdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(0, fdof);
      err = PetscSectionGetFieldDof(section, v, 1, &fdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fdof);
    } else {
      err = PetscSectionGetFieldDof(section, v, 0, &fdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(spaceDim, fdof);
      err = PetscSectionGetFieldDof(section, v, 1, &fdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(0, fdof);
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

  fault->verifyConfiguration(*mesh);
} // _initialize

// ----------------------------------------------------------------------
// Determine if vertex is a Lagrange multiplier constraint vertex.
bool
pylith::faults::TestFaultCohesiveKin::_isConstraintVertex(const int vertex) const
{ // _isConstraintVertex
  CPPUNIT_ASSERT(_data);

  const int numFaultVertices = _data->numFaultVertices;
  bool isFound = false;
  for (int i=0; i < _data->numFaultVertices; ++i)
    if (_data->verticesLagrange[i] == vertex) {
      isFound = true;
      break;
    } // if
  return isFound;
} // _isConstraintVertex


// End of file 
