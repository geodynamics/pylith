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

#include "TestAbsorbingDampers.hh" // Implementation of class methods

#include "pylith/bc/AbsorbingDampers.hh" // USES AbsorbingDampers

#include "data/AbsorbingDampersData.hh" // USES AbsorbingDampersData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES std::runtime_erro

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampers );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestAbsorbingDampers::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestAbsorbingDampers::testConstructor(void)
{ // testConstructor
  AbsorbingDampers bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test db().
void
pylith::bc::TestAbsorbingDampers::testDB(void)
{ // testDB
  const std::string& label = "my db";
  spatialdata::spatialdb::SimpleDB db(label.c_str());
  AbsorbingDampers bc;
  bc.db(&db);
  
  CPPUNIT_ASSERT(0 != bc._db);
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._db->label()));
} // testDB
    
// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestAbsorbingDampers::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  DM             subMesh = boundaryMesh.dmMesh();
  PetscInt       cStart, cEnd, vStart, vEnd;
  PetscErrorCode err;

  // Check boundary mesh
  CPPUNIT_ASSERT(subMesh);

  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(subMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  const int cellDim = boundaryMesh.dimension();
  const int numCorners = _data->numCorners;
  const int spaceDim = _data->spaceDim;
  const int numVertices = vEnd-vStart;
  const int numCells = cEnd-cStart;
  //const int boundaryDepth = submesh->depth()-1; // depth of boundary cells

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  PetscInt dp = 0;
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone;
    PetscInt        numCorners;

    // Assume non-interpolated mesh
    err = DMPlexGetConeSize(subMesh, c, &numCorners);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);
    err = DMPlexGetCone(subMesh, c, &cone);CHECK_PETSC_ERROR(err);
    for(PetscInt p = 0; p < numCorners; ++p, ++dp) {
      CPPUNIT_ASSERT_EQUAL(_data->cells[dp], cone[p]);
    }
  } // for

  // Check damping constants
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  PetscInt index = 0;
  CPPUNIT_ASSERT(0 != bc._parameters);
  PetscSection dampersSection = bc._parameters->get("damping constants").petscSection();
  Vec          dampersVec     = bc._parameters->get("damping constants").localVector();
  CPPUNIT_ASSERT(dampersSection);CPPUNIT_ASSERT(dampersVec);

  const PylithScalar tolerance = 1.0e-06;
  PetscScalar       *dampersValues;
  err = VecGetArray(dampersVec, &dampersValues);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt dof, off;

    err = PetscSectionGetDof(dampersSection, c, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampersSection, c, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    for(PetscInt iQuad=0; iQuad < numQuadPts; ++iQuad)
      for(PetscInt iDim =0; iDim < spaceDim; ++iDim) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dampersValues[off+iQuad*spaceDim+iDim]/_data->dampingConsts[index], tolerance);
        ++index;
      } // for
  } // for
  err = VecRestoreArray(dampersVec, &dampersValues);CHECK_PETSC_ERROR(err);
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestAbsorbingDampers::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  DM             subMesh = boundaryMesh.dmMesh();
  PetscErrorCode err;

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const PylithScalar t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  DM             dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  const PylithScalar* valsE = _data->valsResidual;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const int totalNumVertices = vEnd - vStart;
  const int sizeE = _data->spaceDim * totalNumVertices;

  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  PetscScalar *vals;
  PetscInt     size;

  CPPUNIT_ASSERT(residualSection);CPPUNIT_ASSERT(residualVec);
  err = PetscSectionGetStorageSize(residualSection, &size);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual->view("RESIDUAL");

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
  err = VecRestoreArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  DM subMesh = boundaryMesh.dmMesh();

  topology::Field<topology::Mesh>& solution = fields.solution();
  topology::Jacobian jacobian(solution);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.assemble("final_assembly");

  DM             dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const PylithScalar* valsE = _data->valsJacobian;
  const int totalNumVertices = vEnd - vStart;
  const int nrowsE = totalNumVertices * _data->spaceDim;
  const int ncolsE = totalNumVertices * _data->spaceDim;

  const PetscMat jacobianMat = jacobian.matrix();
  int nrows = 0;
  int ncols = 0;
  err = MatGetSize(jacobianMat, &nrows, &ncols);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  err = MatConvert(jacobianMat, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);CHECK_PETSC_ERROR(err);

  scalar_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  err = MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);CHECK_PETSC_ERROR(err);

#if 0
  std::cout << "JACOBIAN\n";
  for (int iRow=0, i=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol, ++i)
      std::cout << "  iRow: " << iRow << ", iCol: " << iCol << ", value: " << vals[i] << ", valueE: " << valsE[i] << std::endl;
#endif

  const PylithScalar tolerance = 1.0e-06;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valsE[index], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], vals[index], tolerance);
    } // for
  err = MatDestroy(&jDense);CHECK_PETSC_ERROR(err);
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobianLumped().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  DM subMesh = boundaryMesh.dmMesh();
  CPPUNIT_ASSERT(subMesh);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.complete();

  DM             dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const PylithScalar* valsMatrixE = _data->valsJacobian;
  const int totalNumVertices = vEnd - vStart;
  const int sizeE = totalNumVertices * _data->spaceDim;
  scalar_array valsE(sizeE);
  const int spaceDim = _data->spaceDim;
  for (int iVertex=0; iVertex < totalNumVertices; ++iVertex)
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const int indexRow = (iVertex*spaceDim+iDim)*totalNumVertices*spaceDim;
      PylithScalar value = 0.0;
      for (int jVertex=0; jVertex < totalNumVertices; ++jVertex)
	value += valsMatrixE[indexRow + jVertex*spaceDim+iDim];
      valsE[iVertex*spaceDim+iDim] = value;
    } // for

#if 0 // DEBUGGING
  jacobian.view("JACOBIAN");
  std::cout << "\n\nJACOBIAN FULL" << std::endl;
  const int n = totalNumVertices*spaceDim;
  for (int r=0; r < n; ++r) {
    for (int c=0; c < n; ++c) 
      std::cout << "  " << valsMatrixE[r*n+c];
    std::cout << "\n";
  } // for
#endif // DEBUGGING

  PetscSection jacobianSection = jacobian.petscSection();
  Vec          jacobianVec     = jacobian.localVector();
  PetscScalar *vals;
  PetscInt     size;

  CPPUNIT_ASSERT(jacobianSection);CPPUNIT_ASSERT(jacobianVec);
  err = PetscSectionGetStorageSize(jacobianSection, &size);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(jacobianVec, &vals);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
  err = VecRestoreArray(jacobianVec, &vals);CHECK_PETSC_ERROR(err);
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_initialize(topology::Mesh* mesh,
					      AbsorbingDampers* const bc,
					      topology::SolutionFields* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);

  try {
    // Setup mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);

    // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;
    cs.setSpaceDim(mesh->dimension());
    cs.initialize();
    mesh->coordsys(&cs);
    mesh->nondimensionalize(normalizer);

    // Set up quadrature
    _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			    _data->basisDerivRef, _data->numQuadPts, 
			    _data->numBasis, _data->cellDim,
			    _data->quadPts, _data->numQuadPts, _data->cellDim,
			    _data->quadWts, _data->numQuadPts,
			    _data->spaceDim);

    // Set up database
    spatialdata::spatialdb::SimpleDB db("TestAbsorbingDampers");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

    const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

    bc->quadrature(_quadrature);
    bc->timeStep(_data->dt);
    bc->label(_data->label);
    bc->db(&db);
    bc->createSubMesh(*mesh);
    bc->initialize(*mesh, upDir);

    //bc->_boundaryMesh->view("BOUNDARY MESH");

    // Setup fields
    CPPUNIT_ASSERT(0 != fields);
    fields->add("residual", "residual");
    fields->add("dispIncr(t->t+dt)", "displacement_increment");
    fields->add("disp(t)", "displacement");
    fields->add("disp(t-dt)", "displacement");
    fields->add("velocity(t)", "velocity");
    fields->solutionName("dispIncr(t->t+dt)");

    topology::Field<topology::Mesh>& residual = fields->get("residual");
    DM             dmMesh = mesh->dmMesh();
    PetscInt       vStart, vEnd;
    PetscErrorCode err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    residual.newSection(pylith::topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
    residual.allocate();
    residual.zero();
    fields->copyLayout("residual");

    const int totalNumVertices = vEnd - vStart;
    const int fieldSize = _data->spaceDim * totalNumVertices;
    PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
    Vec          dispTIncrVec     = fields->get("dispIncr(t->t+dt)").localVector();
    PetscSection dispTSection     = fields->get("disp(t)").petscSection();
    Vec          dispTVec         = fields->get("disp(t)").localVector();
    PetscSection dispTmdtSection  = fields->get("disp(t-dt)").petscSection();
    Vec          dispTmdtVec      = fields->get("disp(t-dt)").localVector();
    PetscSection velSection       = fields->get("velocity(t)").petscSection();
    Vec          velVec           = fields->get("velocity(t)").localVector();
    PetscScalar *dispTIncrArray, *dispTArray, *dispTmdtArray, *velArray;
    const int spaceDim = _data->spaceDim;
    const PylithScalar dt = _data->dt;

    CPPUNIT_ASSERT(dispTIncrSection);CPPUNIT_ASSERT(dispTIncrVec);
    CPPUNIT_ASSERT(dispTSection);CPPUNIT_ASSERT(dispTVec);
    CPPUNIT_ASSERT(dispTmdtSection);CPPUNIT_ASSERT(dispTmdtVec);
    CPPUNIT_ASSERT(velSection);CPPUNIT_ASSERT(velVec);
    err = VecGetArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(dispTVec,     &dispTArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(dispTmdtVec,  &dispTmdtArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(velVec,       &velArray);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt dof, off;

      err = PetscSectionGetDof(dispTIncrSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTIncrSection, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        dispTIncrArray[off+d] = _data->fieldTIncr[(v - vStart)*spaceDim+d];
      }
      err = PetscSectionGetDof(dispTSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTSection, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        dispTIncrArray[off+d] = _data->fieldT[(v - vStart)*spaceDim+d];
      }
      err = PetscSectionGetDof(dispTmdtSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(dispTmdtSection, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        dispTIncrArray[off+d] = _data->fieldTmdt[(v - vStart)*spaceDim+d];
      }
      err = PetscSectionGetDof(velSection, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(velSection, v, &off);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        velArray[off+d] = (_data->fieldTIncr[(v - vStart)*spaceDim+d] +
                           _data->fieldT[(v - vStart)*spaceDim+d] -
                           _data->fieldTmdt[(v - vStart)*spaceDim+d]) / (2*dt);
      }
    } // for
    err = VecRestoreArray(dispTIncrVec, &dispTIncrArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(dispTVec,     &dispTArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(dispTmdtVec,  &dispTmdtArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(velVec,       &velArray);CHECK_PETSC_ERROR(err);
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
