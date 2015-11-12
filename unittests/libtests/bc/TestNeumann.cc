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

#include "TestNeumann.hh" // Implementation of class methods

#include "pylith/bc/Neumann.hh" // USES Neumann

#include "data/NeumannData.hh" // USES NeumannData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "data/NeumannDataQuad4.hh" // USES NeumannDataQuad4
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

#include <stdexcept> // USES std::runtime_erro

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumann );

// ----------------------------------------------------------------------
namespace pylith {
  namespace bc {
    namespace _TestNeumann {
      const int ncells = 2;
      const int numQuadPts = 2;
      const int spaceDim = 2;

      const PylithScalar initial[ncells*numQuadPts*spaceDim] = {
	0.3,  0.4,    0.7,  0.6,
	1.3,  1.4,    1.7,  1.6,
      };
      const PylithScalar rate[ncells*numQuadPts*spaceDim] = {
	-0.2,  -0.1,   0.4, 0.3,
	-1.2,  -1.1,   1.4, 1.3,
      };
      const PylithScalar rateTime[ncells*numQuadPts] = {
	0.5,   0.8,
	0.6,   0.9,
      };
      const PylithScalar change[ncells*numQuadPts*spaceDim] = {
	1.3,  1.4,    1.7,  1.6,
	2.3,  2.4,    2.7,  2.6,
      };
      const PylithScalar changeTime[ncells*numQuadPts] = {
	2.0,  2.4,
	2.1,  2.5,
      };

      const PylithScalar tValue = 2.2;
      const PylithScalar valuesRate[ncells*numQuadPts*spaceDim] = {
	-0.34,  -0.17,  0.56,   0.42,
	-1.92,  -1.76,  1.82,   1.69,
      };
      const PylithScalar valuesChange[ncells*numQuadPts*spaceDim] = {
	1.3,  1.4,   0.0,  0.0,
	2.3,  2.4,   0.0,  0.0,
      };
      const PylithScalar valuesChangeTH[ncells*numQuadPts*spaceDim] = {
	1.3*0.98,  1.4*0.98,    0.0,  0.0,
	2.3*0.99,  2.4*0.99,    0.0,  0.0,
      };

      // Check values in section against expected values.
      static
      void _checkValues(const PylithScalar* valuesE,
			const int fiberDimE,
			const topology::Field& field);
    } // _TestNeumann
  } // bc
} // pylith

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumann::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _quadrature = new feassemble::Quadrature();CPPUNIT_ASSERT(_quadrature);

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestNeumann::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestNeumann::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  Neumann bc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test _getLabel().
void
pylith::bc::TestNeumann::test_getLabel(void)
{ // test_getLabel
  PYLITH_METHOD_BEGIN;

  Neumann bc;
  
  const std::string& label = "traction bc";
  bc.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._getLabel()));

  PYLITH_METHOD_END;
} // test_getLabel

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestNeumann::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::Mesh& boundaryMesh = *bc._boundaryMesh;

  PetscDM subMesh = boundaryMesh.dmMesh();assert(subMesh);
  topology::Stratum verticesStratum(subMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const int numVertices = verticesStratum.size();

  topology::Stratum cellsStratum(subMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  const PetscInt numCells = cellsStratum.size();

  const int cellDim = boundaryMesh.dimension();
  const int numCorners = _data->numCorners;
  const int spaceDim = _data->spaceDim;

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  PetscInt dp = 0;
  PetscErrorCode err = 0;
  for (PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt  vertices[32];
    PetscInt *closure = NULL;
    PetscInt  closureSize, numCorners = 0;

    err = DMPlexGetTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    for(PetscInt p = 0; p < closureSize*2; p += 2) {
      const PetscInt point = closure[p];
      if ((point >= vStart) && (point < vEnd)) {
        vertices[numCorners++] = point;
      } // if
    } // for
    err = DMPlexRestoreTransitiveClosure(subMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);
    if (boundaryMesh.dimension() > 1) {
      PetscInt first;
      for (first = 0; first < numCorners; ++first) if (_data->cells[dp] == vertices[first]) break;
      CPPUNIT_ASSERT(first < numCorners);
      for (PetscInt p = 0; p < numCorners; ++p, ++dp)
        CPPUNIT_ASSERT_EQUAL(_data->cells[dp], vertices[(p+first)%numCorners]);
    } else {
      for(PetscInt p = 0; p < numCorners; ++p, ++dp) {
        CPPUNIT_ASSERT_EQUAL(_data->cells[dp], vertices[p]);
      } // for
    }
  } // for

  // Check traction values
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  CPPUNIT_ASSERT(bc._parameters);
  topology::VecVisitorMesh initialVisitor(bc._parameters->get("initial"));
  const PetscScalar* initialArray = initialVisitor.localArray();CPPUNIT_ASSERT(initialArray);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar pressureScale = _data->pressureScale;

  for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    const PetscInt off = initialVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(numQuadPts*spaceDim, initialVisitor.sectionDof(c));

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int iDim =0; iDim < spaceDim; ++iDim, ++index) {
        const PylithScalar tractionE = _data->tractionsCell[index];
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionE, initialArray[off+iQuad*spaceDim+iDim]*pressureScale, tolerance);
      } // for
  } // for

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field& residual = fields.get("residual");
  const PylithScalar t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  const int totalNumVertices = verticesStratum.size();

  const PylithScalar* residualE = _data->valsResidual;
  const int spaceDim = _data->spaceDim;

  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();
  //residual.view("RESIDUAL");

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar residualScale = _data->pressureScale * pow(_data->lengthScale, _data->spaceDim-1);

  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = residualVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, residualVisitor.sectionDof(v));

    for (int iDim=0; iDim < spaceDim; ++iDim, ++index) {
      if (fabs(residualE[index]) > 1.0) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+iDim]/residualE[index]*residualScale, tolerance);
      } else {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(residualE[index], residualArray[off+iDim]*residualScale, tolerance);
      } // if/else
    } // for
  } // for

  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test _queryDatabases().
void
pylith::bc::TestNeumann::test_queryDatabases(void)
{ // test_queryDatabases
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const PylithScalar pressureScale = _data->pressureScale;
  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);

  // bc._parameters->view("PARAMETERS"); // DEBUGGING

  // Check initial values.
  const topology::Field& initial = bc._parameters->get("initial");
  _TestNeumann::_checkValues(_TestNeumann::initial, numQuadPts*spaceDim, initial);

  // Check rate values.
  const topology::Field& rate = bc._parameters->get("rate");
  _TestNeumann::_checkValues(_TestNeumann::rate, numQuadPts*spaceDim, rate);

  // Check rate start time.
  const topology::Field& rateTime = bc._parameters->get("rate time");
  _TestNeumann::_checkValues(_TestNeumann::rateTime, numQuadPts, rateTime);

  // Check change values.
  const topology::Field& change = bc._parameters->get("change");
  _TestNeumann::_checkValues(_TestNeumann::change, numQuadPts*spaceDim, change);

  // Check change start time.
  const topology::Field& changeTime = bc._parameters->get("change time");
  _TestNeumann::_checkValues(_TestNeumann::changeTime, numQuadPts, changeTime);
  th.close();

  PYLITH_METHOD_END;
} // test_queryDatabases

// ----------------------------------------------------------------------
// Test _paramsLocalToGlobal().
void
pylith::bc::TestNeumann::test_paramsLocalToGlobal(void)
{ // test_paramsLocalToGlobal
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);

  const PylithScalar pressureScale = _data->pressureScale;
  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  bc._paramsLocalToGlobal(upDir);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);

  // Orientation for quad4 is +x, -y for shear and normal tractions.
  CPPUNIT_ASSERT_EQUAL(2, spaceDim); 
  const int ncells = _TestNeumann::ncells;
  scalar_array valuesE(ncells*numQuadPts*spaceDim);
  
  // Check initial values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::initial[i+0]; // x
    valuesE[i+1] = -_TestNeumann::initial[i+1]; // y
  } // for
  const topology::Field& initial = bc._parameters->get("initial");
  _TestNeumann::_checkValues(&valuesE[0], numQuadPts*spaceDim, initial);

  // Check rate values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::rate[i+0]; // x
    valuesE[i+1] = -_TestNeumann::rate[i+1]; // y
  } // for
  const topology::Field& rate = bc._parameters->get("rate");
  _TestNeumann::_checkValues(&valuesE[0], numQuadPts*spaceDim, rate);

  // Check change values.
  for (int i=0; i < valuesE.size(); i+=spaceDim) {
    valuesE[i+0] = _TestNeumann::change[i+0]; // x
    valuesE[i+1] = -_TestNeumann::change[i+1]; // y
  } // for
  const topology::Field& change = bc._parameters->get("change");
  _TestNeumann::_checkValues(&valuesE[0], numQuadPts*spaceDim, change);

  PYLITH_METHOD_END;
} // test_paramsLocalToGlobal

// ----------------------------------------------------------------------
// Test _calculateValue() with initial value.
void
pylith::bc::TestNeumann::test_calculateValueInitial(void)
{ // test_calculateValueInitial
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbInitial(&dbInitial);

  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);
  
  // Check values.
  const topology::Field& value = bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::initial, numQuadPts*spaceDim, value);

  PYLITH_METHOD_END;
} // test_calculateValueInitial

// ----------------------------------------------------------------------
// Test _calculateValue() with rate.
void
pylith::bc::TestNeumann::test_calculateValueRate(void)
{ // test_calculateValueRate
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbRate(&dbRate);

  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);
  
  // Check values.
  const topology::Field& value = bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesRate, numQuadPts*spaceDim, value);

  PYLITH_METHOD_END;
} // test_calculateValueRate

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change.
void
pylith::bc::TestNeumann::test_calculateValueChange(void)
{ // test_calculateValueChange
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  bc.dbChange(&dbChange);

  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);
  
  // Check values.
  const topology::Field& value = bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesChange, numQuadPts*spaceDim, value);

  PYLITH_METHOD_END;
} // test_calculateValueChange

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change w/time history.
void
pylith::bc::TestNeumann::test_calculateValueChangeTH(void)
{ // test_calculateValueChangeTH
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);
  
  // Check values.
  const topology::Field& value = bc._parameters->get("value");
  _TestNeumann::_checkValues(_TestNeumann::valuesChangeTH, numQuadPts*spaceDim, value);

  PYLITH_METHOD_END;
} // test_calculateValueChangeTH

// ----------------------------------------------------------------------
// Test _calculateValue() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestNeumann::test_calculateValueAll(void)
{ // test_calculateValueAll
  PYLITH_METHOD_BEGIN;

  delete _data; _data = new NeumannDataQuad4();
  feassemble::GeometryLine2D geometry;
  CPPUNIT_ASSERT(_quadrature);
  _quadrature->refGeometry(&geometry);

  topology::Mesh mesh;
  Neumann bc;
  _preinitialize(&mesh, &bc);

  spatialdata::spatialdb::SimpleDB dbInitial("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/quad4_traction_initial.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/quad4_traction_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("_TestNeumann _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/quad4_traction_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("_TestNeumann _queryDatabases");
  th.filename("data/quad4_traction.timedb");

  bc.dbInitial(&dbInitial);
  bc.dbRate(&dbRate);
  bc.dbChange(&dbChange);
  bc.dbTimeHistory(&th);

  const PylithScalar timeScale = _data->timeScale;
  bc._queryDatabases();
  bc._calculateValue(_TestNeumann::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int spaceDim = _TestNeumann::spaceDim;
  const int numQuadPts = _TestNeumann::numQuadPts;
  CPPUNIT_ASSERT(bc._parameters);
  
  // Check values.
  const int ncells = _TestNeumann::ncells;
  scalar_array valuesE(ncells*numQuadPts*spaceDim);
  for (int i=0; i < valuesE.size(); ++i)
    valuesE[i] = 
      _TestNeumann::initial[i] +
      _TestNeumann::valuesRate[i] +
      _TestNeumann::valuesChangeTH[i];
  
  const topology::Field& value = bc._parameters->get("value");
  _TestNeumann::_checkValues(&valuesE[0], numQuadPts*spaceDim, value);

  PYLITH_METHOD_END;
} // test_calculateValueAll

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_preinitialize(topology::Mesh* mesh,
					Neumann* const bc) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(bc);

  // Set up mesh
  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  // Set up coordinates
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

  bc->quadrature(_quadrature);
  bc->label(_data->label);
  bc->normalizer(normalizer);
  bc->createSubMesh(*mesh);

  PYLITH_METHOD_END;
} // _preinitialize

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_initialize(topology::Mesh* mesh,
				     Neumann* const bc,
				     topology::SolutionFields* fields) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(bc);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_quadrature);

  _preinitialize(mesh, bc);

  // Set up database
  spatialdata::spatialdb::SimpleDB db("TestNeumann");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->spatialDBFilename);
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

  bc->dbInitial(&db);
  bc->initialize(*mesh, upDir);

  // Set up fields
  CPPUNIT_ASSERT(fields);
  fields->add("residual", "residual");
  fields->add("disp(t), bc(t+dt)", "displacement");
  fields->solutionName("disp(t), bc(t+dt)");

  topology::Field& residual = fields->get("residual");
  residual.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  residual.allocate();
  residual.scale(_data->lengthScale);
  residual.zero();

  fields->copyLayout("residual");

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Check values in section against expected values.
void
pylith::bc::_TestNeumann::_checkValues(const PylithScalar* valuesE,
				       const int fiberDimE,
				       const topology::Field& field)
{ // _checkValues
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(valuesE);

  const topology::Mesh& boundaryMesh = field.mesh();
  const PetscInt ncells = _TestNeumann::ncells;

  PetscDM subMesh = boundaryMesh.dmMesh();CPPUNIT_ASSERT(subMesh);
  topology::Stratum cellsStratum(subMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  CPPUNIT_ASSERT_EQUAL(ncells, cellsStratum.size());

  const PylithScalar scale = field.scale();

  // Check values associated with BC.
  topology::VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();

  const PylithScalar tolerance = 1.0e-06;
  for (PetscInt c = cStart, icell = 0; c < cEnd; ++c, ++icell) {
    const PetscInt off = fieldVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(c));

    for (int iDim=0; iDim < fiberDimE; ++iDim) {
      if (valuesE[icell*fiberDimE+iDim] != 0.0) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldArray[off+iDim]/valuesE[icell*fiberDimE+iDim]*scale, tolerance);
      } else {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[icell*fiberDimE+iDim], fieldArray[off+iDim]*scale, tolerance);
      } // if/else
    } // for
  } // for

  PYLITH_METHOD_END;
} // _checkValues


// End of file 
