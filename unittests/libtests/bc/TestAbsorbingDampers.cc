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
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(_quadrature);
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
  
  CPPUNIT_ASSERT(bc._db);
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

  CPPUNIT_ASSERT(_data);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM subMesh = boundaryMesh.dmMesh();
  PetscInt cStart, cEnd, vStart, vEnd;
  PetscErrorCode err = 0;

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
    PetscInt *closure = NULL;
    PetscInt closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        closure[numCorners++] = point;
      }
    }
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);
    for(PetscInt p = 0; p < numCorners; ++p, ++dp) {
      CPPUNIT_ASSERT_EQUAL(_data->cells[dp], closure[p]);
    } // for
    err = DMPlexRestoreTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);CHECK_PETSC_ERROR(err);
  } // for

  // Check damping constants
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  PetscInt index = 0;
  CPPUNIT_ASSERT(bc._parameters);
  PetscSection dampingConstsSection = bc._parameters->get("damping constants").petscSection();
  PetscVec dampingConstsVec = bc._parameters->get("damping constants").localVector();
  CPPUNIT_ASSERT(dampingConstsSection);CPPUNIT_ASSERT(dampingConstsVec);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar* dampingConstsE = _data->dampingConsts;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  PetscScalar* dampingConstsArray;
  err = VecGetArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt dof, off;

    err = PetscSectionGetDof(dampingConstsSection, c, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dampingConstsSection, c, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    for(PetscInt iQuad=0; iQuad < numQuadPts; ++iQuad)
      for(PetscInt iDim =0; iDim < spaceDim; ++iDim) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dampingConstsArray[off+iQuad*spaceDim+iDim]/dampingConstsE[index]*dampingConstsScale, tolerance);
        ++index;
      } // for
  } // for
  err = VecRestoreArray(dampingConstsVec, &dampingConstsArray);CHECK_PETSC_ERROR(err);
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestAbsorbingDampers::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM             subMesh = boundaryMesh.dmMesh();
  PetscErrorCode err;

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const PylithScalar t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  PetscDM dmMesh = mesh.dmMesh();
  PetscInt vStart, vEnd;
  const PylithScalar* valsE = _data->valsResidual;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar velocityScale = 1.0; // Input velocity is nondimensional.
  const PylithScalar residualScale = dampingConstsScale*velocityScale*pow(_data->lengthScale, _data->spaceDim-1);

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const int totalNumVertices = vEnd - vStart;
  const int sizeE = _data->spaceDim * totalNumVertices;

  PetscSection residualSection = residual.petscSection();
  PetscVec residualVec = residual.localVector();
  PetscScalar *vals;
  PetscInt size;

  CPPUNIT_ASSERT(residualSection);CPPUNIT_ASSERT(residualVec);
  err = PetscSectionGetStorageSize(residualSection, &size);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual->view("RESIDUAL");

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i]*residualScale, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i]/residualScale, vals[i], tolerance);
  err = VecRestoreArray(residualVec, &vals);CHECK_PETSC_ERROR(err);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM subMesh = boundaryMesh.dmMesh();

  topology::Field<topology::Mesh>& solution = fields.solution();
  topology::Jacobian jacobian(solution);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.assemble("final_assembly");

  PetscDM dmMesh = mesh.dmMesh();
  PetscInt vStart, vEnd;
  PetscErrorCode err = 0;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const PylithScalar* valsE = _data->valsJacobian;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar jacobianScale = dampingConstsScale*pow(_data->lengthScale, _data->spaceDim-1);

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

#if 0 // DEBUGGING
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
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valsE[index]*jacobianScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index]/jacobianScale, vals[index], tolerance);
    } // for
  err = MatDestroy(&jDense);CHECK_PETSC_ERROR(err);
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobianLumped().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  CPPUNIT_ASSERT(_data);

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
  PetscDM subMesh = boundaryMesh.dmMesh();
  CPPUNIT_ASSERT(subMesh);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.complete();

  PetscDM dmMesh = mesh.dmMesh();
  PetscInt vStart, vEnd;
  PetscErrorCode err = 0;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const PylithScalar* valsMatrixE = _data->valsJacobian;
  const int totalNumVertices = vEnd - vStart;
  const int sizeE = totalNumVertices * _data->spaceDim;
  scalar_array valsE(sizeE);
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar jacobianScale = dampingConstsScale*pow(_data->lengthScale, _data->spaceDim-1);

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
  PetscVec jacobianVec = jacobian.localVector();
  PetscScalar *vals;
  PetscInt size;

  CPPUNIT_ASSERT(jacobianSection);CPPUNIT_ASSERT(jacobianVec);
  err = PetscSectionGetStorageSize(jacobianSection, &size);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(jacobianVec, &vals);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i]*jacobianScale, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i]/jacobianScale, vals[i], tolerance);
  err = VecRestoreArray(jacobianVec, &vals);CHECK_PETSC_ERROR(err);
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_initialize(topology::Mesh* mesh,
					      AbsorbingDampers* const bc,
					      topology::SolutionFields* fields) const
{ // _initialize
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(bc);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  try {
    // Setup mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);

    // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(_data->lengthScale);
    normalizer.pressureScale(_data->pressureScale);
    normalizer.densityScale(_data->densityScale);
    normalizer.timeScale(_data->timeScale);
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
    bc->normalizer(normalizer);
    bc->createSubMesh(*mesh);
    bc->initialize(*mesh, upDir);

    //bc->_boundaryMesh->view("BOUNDARY MESH");

    // Setup fields
    CPPUNIT_ASSERT(fields);
    fields->add("residual", "residual");
    fields->add("dispIncr(t->t+dt)", "displacement_increment");
    fields->add("disp(t)", "displacement");
    fields->add("disp(t-dt)", "displacement");
    fields->add("velocity(t)", "velocity");
    fields->solutionName("dispIncr(t->t+dt)");

    topology::Field<topology::Mesh>& residual = fields->get("residual");
    PetscDM             dmMesh = mesh->dmMesh();
    PetscInt       vStart, vEnd;
    PetscErrorCode err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    residual.newSection(pylith::topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
    residual.allocate();
    residual.zero();
    residual.scale(normalizer.lengthScale());
    fields->copyLayout("residual");

    const int totalNumVertices = vEnd - vStart;
    const int fieldSize = _data->spaceDim * totalNumVertices;
    PetscSection dispTIncrSection = fields->get("dispIncr(t->t+dt)").petscSection();
    PetscVec dispTIncrVec = fields->get("dispIncr(t->t+dt)").localVector();
    PetscSection dispTSection = fields->get("disp(t)").petscSection();
    PetscVec dispTVec = fields->get("disp(t)").localVector();
    PetscSection dispTmdtSection  = fields->get("disp(t-dt)").petscSection();
    PetscVec dispTmdtVec = fields->get("disp(t-dt)").localVector();
    PetscSection velSection = fields->get("velocity(t)").petscSection();
    PetscVec velVec = fields->get("velocity(t)").localVector();
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
