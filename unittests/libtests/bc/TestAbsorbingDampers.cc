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
// Copyright (c) 2010-2014 University of California, Davis
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
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampers );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _quadrature = new feassemble::Quadrature();CPPUNIT_ASSERT(_quadrature);

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestAbsorbingDampers::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestAbsorbingDampers::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  AbsorbingDampers bc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test db().
void
pylith::bc::TestAbsorbingDampers::testDB(void)
{ // testDB
  PYLITH_METHOD_BEGIN;

  const std::string& label = "my db";
  spatialdata::spatialdb::SimpleDB db(label.c_str());
  AbsorbingDampers bc;
  bc.db(&db);
  
  CPPUNIT_ASSERT(bc._db);
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._db->label()));

  PYLITH_METHOD_END;
} // testDB
    
// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestAbsorbingDampers::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(_data);

  const topology::Mesh& boundaryMesh = *bc._boundaryMesh;

  PetscDM subMesh = boundaryMesh.dmMesh();CPPUNIT_ASSERT(subMesh);

  // Cells
  topology::Stratum cellsStratum(subMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  // Vertices
  topology::Stratum verticesStratum(subMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const PetscInt numVertices = verticesStratum.size();

  const int cellDim = boundaryMesh.dimension();
  const int numCorners = _data->numCorners;
  const int spaceDim = _data->spaceDim;

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  PetscErrorCode err = 0;
  PetscInt dp = 0;
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt  vertices[32];
    PetscInt *closure = NULL;
    PetscInt closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        vertices[numCorners++] = point;
      }
    }
    err = DMPlexRestoreTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);
    // Allow cyclic permutations to handle both interpolated and non-interpolated cases
    if (boundaryMesh.dimension() > 1) {
      PetscInt first;
      for (first = 0; first < numCorners; ++first) if (_data->cells[dp] == vertices[first]) break;
      CPPUNIT_ASSERT(first < numCorners);
      for (PetscInt p = 0; p < numCorners; ++p, ++dp)
        CPPUNIT_ASSERT_EQUAL(_data->cells[dp], vertices[(p+first)%numCorners]);
    } else {
      for(PetscInt p = 0; p < numCorners; ++p, ++dp) {
        CPPUNIT_ASSERT_EQUAL(_data->cells[dp], vertices[p]);
      }
    } // for
  } // for

  // Check damping constants
  CPPUNIT_ASSERT(bc._parameters);
  topology::VecVisitorMesh dampingConstsVisitor(bc._parameters->get("damping constants"));
  const PetscScalar* dampingConstsArray = dampingConstsVisitor.localArray();

  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  const PylithScalar* dampingConstsE = _data->dampingConsts;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt c = cStart, index=0; c < cEnd; ++c) {
    const PetscInt off = dampingConstsVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dampingConstsVisitor.sectionDof(c));
    for(PetscInt iQuad=0; iQuad < numQuadPts; ++iQuad)
      for(PetscInt iDim =0; iDim < spaceDim; ++iDim, ++index) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dampingConstsArray[off+iQuad*spaceDim+iDim]/dampingConstsE[index]*dampingConstsScale, tolerance);
      } // for
  } // for

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestAbsorbingDampers::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::Mesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM             subMesh = boundaryMesh.dmMesh();
  PetscErrorCode err;

  topology::Field& residual = fields.get("residual");
  const PylithScalar t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  PetscDM dmMesh = mesh.dmMesh();
  PetscInt vStart, vEnd;
  const PylithScalar* valsE = _data->valsResidual;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar velocityScale = 1.0; // Input velocity is nondimensional.
  const PylithScalar residualScale = dampingConstsScale*velocityScale*pow(_data->lengthScale, _data->spaceDim-1);

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  const int totalNumVertices = vEnd - vStart;
  const int sizeE = _data->spaceDim * totalNumVertices;

  PetscSection residualSection = residual.localSection();
  PetscVec residualVec = residual.localVector();
  PetscScalar *vals;
  PetscInt size;

  CPPUNIT_ASSERT(residualSection);CPPUNIT_ASSERT(residualVec);
  err = PetscSectionGetStorageSize(residualSection, &size);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual->view("RESIDUAL");

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(residualVec, &vals);PYLITH_CHECK_ERROR(err);
  for(int i = 0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i]*residualScale, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i]/residualScale, vals[i], tolerance);
  err = VecRestoreArray(residualVec, &vals);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::Mesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM subMesh = boundaryMesh.dmMesh();CPPUNIT_ASSERT(subMesh);

  topology::Field& solution = fields.solution();
  topology::Jacobian jacobian(solution);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.assemble("final_assembly");

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const PetscInt totalNumVertices = verticesStratum.size();

  const PylithScalar* valsE = _data->valsJacobian;
  const PylithScalar dampingConstsScale = _data->densityScale * _data->lengthScale / _data->timeScale;
  const PylithScalar jacobianScale = dampingConstsScale*pow(_data->lengthScale, _data->spaceDim-1);

  const int nrowsE = totalNumVertices * _data->spaceDim;
  const int ncolsE = totalNumVertices * _data->spaceDim;
  const PetscMat jacobianMat = jacobian.matrix();CPPUNIT_ASSERT(jacobianMat);
  int nrows = 0;
  int ncols = 0;
  PetscErrorCode err = 0;
  err = MatGetSize(jacobianMat, &nrows, &ncols);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  err = MatConvert(jacobianMat, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);PYLITH_CHECK_ERROR(err);CPPUNIT_ASSERT(jDense);

  scalar_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  err = MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);PYLITH_CHECK_ERROR(err);

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
  err = MatDestroy(&jDense);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobianLumped().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const topology::Mesh& boundaryMesh = *bc._boundaryMesh;
  PetscDM subMesh = boundaryMesh.dmMesh();CPPUNIT_ASSERT(subMesh);

  const PylithScalar t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.complete();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const PetscInt totalNumVertices = verticesStratum.size();

  const PylithScalar* valsMatrixE = _data->valsJacobian;
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

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();
  const PylithScalar tolerance = 1.0e-06;
  for(int v = vStart, index=0; v < vEnd; ++v) {
    const PetscInt off = jacobianVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, jacobianVisitor.sectionDof(v));
    for (int iDim=0; iDim < spaceDim; ++iDim, ++index) {
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobianArray[off+iDim]/valsE[index]*jacobianScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index]/jacobianScale, jacobianArray[off+iDim], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_initialize(topology::Mesh* mesh,
					      AbsorbingDampers* const bc,
					      topology::SolutionFields* fields) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(bc);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  // Setup mesh
  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

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

  topology::Field& residual = fields->get("residual");
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscInt vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  residual.newSection(pylith::topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  residual.allocate();
  residual.zero();
  residual.scale(normalizer.lengthScale());
  fields->copyLayout("residual");

  const int totalNumVertices = vEnd - vStart;
  const int fieldSize = _data->spaceDim * totalNumVertices;

  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::VecVisitorMesh dispTmdtVisitor(fields->get("disp(t-dt)"));
  PetscScalar* dispTmdtArray = dispTmdtVisitor.localArray();

  topology::VecVisitorMesh velocityVisitor(fields->get("velocity(t)"));
  PetscScalar* velocityArray = velocityVisitor.localArray();

  const int spaceDim = _data->spaceDim;
  const PylithScalar dt = _data->dt;

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt doff = v-vStart;

    PetscInt off = dispTIncrVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTIncrVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTIncrArray[off+d] = _data->fieldTIncr[doff*spaceDim+d];
    } // for

    off = dispTVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTArray[off+d] = _data->fieldT[doff*spaceDim+d];
    } // for

    off = dispTmdtVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTmdtVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTmdtArray[off+d] = _data->fieldTmdt[doff*spaceDim+d];
    } // for

    off = velocityVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, velocityVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      velocityArray[off+d] = (_data->fieldTIncr[doff*spaceDim+d] +
			      _data->fieldT[doff*spaceDim+d] -
			      _data->fieldTmdt[doff*spaceDim+d]) / (2*dt);
    } // for
  } // for

  PYLITH_METHOD_END;
} // _initialize


// End of file 
