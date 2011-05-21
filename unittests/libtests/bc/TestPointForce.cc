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
// Copyright (c) 2010 University of California, Davis
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
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/PackedFields.hh" // USES PackedFields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForce );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForce::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestPointForce::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestPointForce::testConstructor(void)
{ // testConstructor
  PointForce bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test normalizer.
void
pylith::bc::TestPointForce::testNormalizer(void)
{ // testNormalizer
  PointForce bc;

  spatialdata::units::Nondimensional normalizer;
  const double scale = 4.0;
  normalizer.lengthScale(4.0);

  bc.normalizer(normalizer);
  CPPUNIT_ASSERT_EQUAL(scale, bc._getNormalizer().lengthScale());
} // testNormalizer

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestPointForce::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);
  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  
  const int numCells = sieveMesh->heightStratum(0)->size();
  const int numForceDOF = _data->numForceDOF;
  const size_t numPoints = _data->numForcePts;

  // Check points
  const int offset = numCells;
  if (numForceDOF > 0) {
    CPPUNIT_ASSERT_EQUAL(numPoints, bc._points.size());
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(_data->forcePoints[i]+offset, bc._points[i]);
  } // if

  CPPUNIT_ASSERT(0 != bc._parameters);
  const ALE::Obj<RealUniformSection>& parametersSection =
    bc._parameters->section();
  CPPUNIT_ASSERT(!parametersSection.isNull());
  const int parametersFiberDim = bc._parameters->fiberDim();
  
  const double tolerance = 1.0e-06;

  // Check values
  const int initialIndex = bc._parameters->sectionIndex("initial");
  const int initialFiberDim = bc._parameters->sectionFiberDim("initial");
  CPPUNIT_ASSERT_EQUAL(numForceDOF, initialFiberDim);

  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;

    CPPUNIT_ASSERT_EQUAL(parametersFiberDim, 
			 parametersSection->getFiberDimension(p_force));
    const double* parametersVertex = parametersSection->restrictPoint(p_force);
    CPPUNIT_ASSERT(parametersVertex);

    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceInitial[i*numForceDOF+iDOF],
				   parametersVertex[initialIndex+iDOF],
				   tolerance);
  } // for

  // Check rate of change
  const int rateIndex = bc._parameters->sectionIndex("rate");
  const int rateFiberDim = bc._parameters->sectionFiberDim("rate");
  CPPUNIT_ASSERT_EQUAL(numForceDOF, rateFiberDim);

  for (int i=0; i < numPoints; ++i) {
    const int p_force = _data->forcePoints[i]+offset;

    CPPUNIT_ASSERT_EQUAL(parametersFiberDim, 
			 parametersSection->getFiberDimension(p_force));
    const double* parametersVertex = parametersSection->restrictPoint(p_force);
    CPPUNIT_ASSERT(parametersVertex);

    for (int iDOF=0; iDOF < numForceDOF; ++iDOF) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->forceRate,
				   parametersVertex[rateIndex+iDOF],
				   tolerance);
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestPointForce::testIntegrateResidual(void)
{ // testIntegrateResidual
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  topology::Field<topology::Mesh> residual(mesh);
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  CPPUNIT_ASSERT(0 != cs);
  const int spaceDim = cs->spaceDim();
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  residual.zero();

  topology::SolutionFields fields(mesh);

  const double t = _data->tResidual;
  bc.integrateResidual(residual, t, &fields);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT(!sieveMesh->depthStratum(0).isNull());

  const double* valsE = _data->residual;
  const int totalNumVertices = sieveMesh->depthStratum(0)->size();
  const int sizeE = spaceDim * totalNumVertices;

  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());
  const double* vals = residualSection->restrictSpace();
  const int size = residualSection->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual.view("RESIDUAL");

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::bc::TestPointForce::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  topology::Mesh mesh;
  PointForce bc;
  _initialize(&mesh, &bc);

  bc.verifyConfiguration(mesh);
} // testVerifyConfiguration

// ----------------------------------------------------------------------
void
pylith::bc::TestPointForce::_initialize(topology::Mesh* mesh,
					 PointForce* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != bc);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

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

  const double upDir[] = { 0.0, 0.0, 1.0 };

  bc->label(_data->label);
  bc->dbInitial(&dbInitial);
  bc->dbRate(&dbRate);
  bc->bcDOF(_data->forceDOF, _data->numForceDOF);
  bc->initialize(*mesh, upDir);
} // _initialize


// End of file 
