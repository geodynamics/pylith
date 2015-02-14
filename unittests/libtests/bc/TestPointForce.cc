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

#include "TestPointForce.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "data/PointForceData.hh" // USES PointForceData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForce );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForce::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestPointForce::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestPointForce::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  PointForce bc;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test normalizer.
void
pylith::bc::TestPointForce::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;

  PointForce bc;

  spatialdata::units::Nondimensional normalizer;
  const double scale = 4.0;
  normalizer.lengthScale(4.0);

  bc.normalizer(normalizer);
  CPPUNIT_ASSERT_EQUAL(scale, bc._getNormalizer().lengthScale());

  PYLITH_METHOD_END;
} // testNormalizer

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestPointForce::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt numCells = cellsStratum.size();

  const int numForceDOF = _data->numForceDOF;
  const size_t numPoints = _data->numForcePts;

  // Check points
  const int offset = numCells;
  if (numForceDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i) {
      CPPUNIT_ASSERT_EQUAL(_data->forcePoints[i]+offset, bc._points[i]);
    } // for
  } // if

  // Check values
  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar forceScale = _data->pressureScale*pow(_data->lengthScale, 2);
  const PylithScalar timeScale = _data->timeScale;

  CPPUNIT_ASSERT(bc._parameters);
  topology::VecVisitorMesh initialVisitor(bc._parameters->get("initial"));
  const PetscScalar* initialArray = initialVisitor.localArray();

  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;

    const PetscInt off = initialVisitor.sectionOffset(p_force);
    CPPUNIT_ASSERT_EQUAL(numForceDOF, initialVisitor.sectionDof(p_force));

    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceInitial[i*numForceDOF+iDOF], initialArray[off+iDOF]*forceScale, tolerance);
    } // for
  } // for

  // Check rate of change
  topology::VecVisitorMesh rateVisitor(bc._parameters->get("rate"));
  const PetscScalar* rateArray = rateVisitor.localArray();

  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;

    const PetscInt off = rateVisitor.sectionOffset(p_force);
    CPPUNIT_ASSERT_EQUAL(numForceDOF, rateVisitor.sectionDof(p_force));
    
    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceRate, rateArray[off+iDOF]*forceScale/timeScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestPointForce::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  topology::Field residual(mesh);
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();CPPUNIT_ASSERT(cs);
  const int spaceDim = cs->spaceDim();
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  residual.zero();

  topology::SolutionFields fields(mesh);

  const PylithScalar t = _data->tResidual/_data->timeScale;
  bc.integrateResidual(residual, t, &fields);

  // residual.view("RESIDUAL"); // DEBUGGING

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const PylithScalar* residualE = _data->residual;

  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar forceScale = _data->pressureScale*pow(_data->lengthScale, 2);
  const PylithScalar residualScale = forceScale;
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
// Test verifyConfiguration().
void
pylith::bc::TestPointForce::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  bc.verifyConfiguration(mesh);

  PYLITH_METHOD_END;
} // testVerifyConfiguration

// ----------------------------------------------------------------------
void
pylith::bc::TestPointForce::_initialize(topology::Mesh* mesh,
					 PointForce* const bc) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(bc);

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

  spatialdata::spatialdb::SimpleDB dbInitial("TestPointForce initial");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename(_data->dbFilename);
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::UniformDB dbRate("TestPointForce rate");
  { // rate db
    const int numValues = 4;
    const char* names[numValues] = {
      "force-rate-x",
      "force-rate-y",
      "force-rate-z",
      "rate-start-time",
    };
    const char* units[numValues] = {
      "newton",
      "newton",
      "newton",
      "s",
    };
    const double values[numValues] = { _data->forceRate,
				       _data->forceRate,
				       _data->forceRate,
				       _data->tRef};
    dbRate.setData(names, units, values, numValues);
  } // rate db

  const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->dbInitial(&dbInitial);
  bc->dbRate(&dbRate);
  bc->bcDOF(_data->forceDOF, _data->numForceDOF);
  bc->normalizer(normalizer);
  bc->initialize(*mesh, upDir);

  PYLITH_METHOD_END;
} // _initialize


// End of file 
