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

#include "TestElasticityImplicitLgDeform.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityImplicitLgDeform.hh" // USES ElasticityImplicitLgDeform
#include "data/IntegratorData.hh" // USES IntegratorData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <math.h> // USES fabs()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _quadrature = new Quadrature();CPPUNIT_ASSERT(_quadrature);
  _data = 0;
  _material = 0;
  _gravityField = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _material; _material = 0;
  delete _gravityField; _gravityField = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestElasticityImplicitLgDeform::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  ElasticityImplicitLgDeform integrator;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void 
pylith::feassemble::TestElasticityImplicitLgDeform::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityImplicitLgDeform integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityImplicitLgDeform::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityImplicitLgDeform integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  topology::Field& residual = fields.get("residual");
  const PylithScalar t = 1.0;
  integrator.integrateResidual(residual, t, &fields);

  const PylithScalar* valsE = _data->valsResidual;

#if 0 // DEBUGGING
  residual.view("RESIDUAL");
  std::cout << "EXPECTED RESIDUAL" << std::endl;
  const int size = _data->numVertices * _data->spaceDim;
  for (int i=0; i < size; ++i)
    std::cout << "  " << valsE[i] << std::endl;
#endif

  const PetscDM dmMesh = mesh.dmMesh();
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, verticesStratum.size());

  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

  const PylithScalar accScale = _data->lengthScale / pow(_data->timeScale, 2);
  const PylithScalar residualScale = _data->densityScale * accScale*pow(_data->lengthScale, _data->spaceDim);

  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = residualVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(_data->spaceDim, residualVisitor.sectionDof(v));

    for (int d=0; d < _data->spaceDim; ++d, ++index) {
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valsE[index]*residualScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], residualArray[off+d]*residualScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityImplicitLgDeform::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityImplicitLgDeform integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);
  integrator._needNewJacobian = true;

  topology::Jacobian jacobian(fields.solution());

  const PylithScalar t = 1.0;
  integrator.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());
  jacobian.assemble("final_assembly");

  const PylithScalar* valsE = _data->valsJacobian;
  const int nrowsE = _data->numVertices * _data->spaceDim;
  const int ncolsE = _data->numVertices * _data->spaceDim;

  const PetscMat jacobianMat = jacobian.matrix();

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

  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 2.0e-05;
  const PylithScalar jacobianScale = _data->densityScale / pow(_data->timeScale, 2) * pow(_data->lengthScale, _data->spaceDim);

  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      const PylithScalar valE = valsE[index];
      if (fabs(valE) > 1.0) {
	// Adjust tolerance based on magnitude of expected value compared to typical Jacobian values of 1.0e+11
	const PylithScalar toleranceAdj = (fabs(valE) < 1.0e+10) ? tolerance*1.0e+11/fabs(valE) : tolerance;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valE*jacobianScale, toleranceAdj);
      } else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, vals[index]*jacobianScale, tolerance);
    } // for
  MatDestroy(&jDense);

  PYLITH_METHOD_END;
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test updateStateVars().
void 
pylith::feassemble::TestElasticityImplicitLgDeform::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityImplicitLgDeform integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const PylithScalar t = 1.0;
  integrator.updateStateVars(t, &fields);

  PYLITH_METHOD_END;
} // testUpdateStateVars

// Initialize elasticity integrator.
void
pylith::feassemble::TestElasticityImplicitLgDeform::_initialize(topology::Mesh* mesh,
								ElasticityImplicitLgDeform* const integrator,
								topology::SolutionFields* fields)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(integrator);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_material);

  const int spaceDim = _data->spaceDim;

  // Setup mesh
  PetscDM dmMesh;

  // Cells and vertices
  const PetscBool interpolate = PETSC_TRUE;
  PetscErrorCode err;
  err = DMPlexCreateFromCellList(PETSC_COMM_WORLD, _data->cellDim, _data->numCells, _data->numVertices, _data->numBasis, interpolate, _data->cells, _data->spaceDim, _data->vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(dmMesh, "domain");

  // Material ids
  PetscInt cStart, cEnd;
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMSetLabelValue(dmMesh, "material-id", c, _data->matId);PYLITH_CHECK_ERROR(err);
  } // for

  // Setup quadrature
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDerivRef, _data->numQuadPts,
			  _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);

  // Setup coordinate system.
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  // Setup scales.
  const PylithScalar timeScale = _data->timeScale;
  const PylithScalar lengthScale = _data->lengthScale;
  const PylithScalar velScale = lengthScale / timeScale;
  const PylithScalar accScale = lengthScale / (timeScale*timeScale);
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Setup material
  spatialdata::spatialdb::SimpleIOAscii iohandler;
  iohandler.filename(_data->matDBFilename);
  spatialdata::spatialdb::SimpleDB dbProperties;
  dbProperties.ioHandler(&iohandler);
  
  _material->id(_data->matId);
  _material->label(_data->matLabel);
  _material->dbProperties(&dbProperties);
  _material->normalizer(normalizer);

  integrator->quadrature(_quadrature);
  integrator->gravityField(_gravityField);
  integrator->timeStep(_data->dt / _data->timeScale);
  integrator->material(_material);
  integrator->initialize(*mesh);

  // Setup fields
  CPPUNIT_ASSERT(fields);
  fields->add("residual", "residual");
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->solutionName("dispIncr(t->t+dt)");
  
  topology::Field& residual = fields->get("residual");
  residual.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR, lengthScale);
  residual.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);

  residual.subfieldsSetup();
  residual.setupSolnChart();
  residual.setupSolnDof(spaceDim);
  residual.allocate();
  residual.zeroAll();
  fields->copyLayout("residual");

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  PetscScalar* dispTArray = dispTVisitor.localArray();CPPUNIT_ASSERT(dispTArray);

  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();CPPUNIT_ASSERT(dispTIncrArray);

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt dtoff = dispTVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTVisitor.sectionDof(v));

    const PetscInt dioff = dispTIncrVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTIncrVisitor.sectionDof(v));

    for(int iDim=0; iDim < spaceDim; ++iDim) {
      dispTArray[dtoff+iDim] = _data->fieldT[iVertex*spaceDim+iDim] / lengthScale;
      dispTIncrArray[dioff+iDim] = _data->fieldTIncr[iVertex*spaceDim+iDim] / lengthScale;
    } // for
  } // for

  PYLITH_METHOD_END;
} // _initialize


// End of file 
