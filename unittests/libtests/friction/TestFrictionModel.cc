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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFrictionModel.hh" // Implementation of class methods

#include "data/StaticFrictionData.hh" // USES StaticFrictionData
#include "data/SlipWeakeningData.hh" // USES SlipWeakeningData

#include "pylith/friction/StaticFriction.hh" // USES StaticFriction
#include "pylith/friction/SlipWeakening.hh" // USES SlipWeakening

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestFrictionModel );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;

// ----------------------------------------------------------------------
// Test label()
void
pylith::friction::TestFrictionModel::testLabel(void)
{ // testLabel
  const std::string& label = "the database";
  StaticFriction friction;
  friction.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction.label()));
} // testLabel
    
// ----------------------------------------------------------------------
// Test timestep()
void
pylith::friction::TestFrictionModel::testTimeStep(void) 
{ // testTimeStep
  const double dt = 2.0;
  StaticFriction friction;
  friction.timeStep(dt);
  
  CPPUNIT_ASSERT_EQUAL(dt, friction.timeStep());
} // testTimeStep

// ----------------------------------------------------------------------
// Test dbProperties()
void
pylith::friction::TestFrictionModel::testDBProperties(void)
{ // testDBProperties
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  StaticFriction friction;
  friction.dbProperties(&db);
  
  CPPUNIT_ASSERT(0 != friction._dbProperties);
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction._dbProperties->label()));
} // testDBProperties

// ----------------------------------------------------------------------
// Test dbStateVars()
void
pylith::friction::TestFrictionModel::testDBStateVars(void)
{ // testDBStateVars
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  StaticFriction friction;
  friction.dbInitialState(&db);
  
  CPPUNIT_ASSERT(0 != friction._dbInitialState);
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction._dbInitialState->label()));
} // testDBStateVars

// ----------------------------------------------------------------------
// Test normalizer()
void
pylith::friction::TestFrictionModel::testNormalizer(void)
{ // testNormalizer
  spatialdata::units::Nondimensional normalizer;
  const double lengthScale = 2.0;
  normalizer.lengthScale(lengthScale);

  StaticFriction friction;
  friction.normalizer(normalizer);
  
  CPPUNIT_ASSERT(0 != friction._normalizer);
  CPPUNIT_ASSERT_EQUAL(lengthScale, friction._normalizer->lengthScale());
} // testNormalizer

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::friction::TestFrictionModel::testInitialize(void)
{ // testInitialize
  const double propertiesE[] = { 0.6, 1000000, 0.4, 1000000 };
  const int numProperties = 2;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = fault.faultMesh().sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const double tolerance = 1.0e-06;

  // Test fieldsPropsStateVars with mesh
  const int fieldsFiberDim = numProperties;
  int index = 0;
  CPPUNIT_ASSERT(0 != friction._fieldsPropsStateVars);
  const ALE::Obj<SubRealUniformSection>& fieldsSection =
    friction._fieldsPropsStateVars->section();
  CPPUNIT_ASSERT(!fieldsSection.isNull());
  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fieldsFiberDim, 
			 fieldsSection->getFiberDimension(*v_iter));
    const double* fieldsVertex = fieldsSection->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(fieldsVertex);
    for (int i = 0; i < numProperties; ++i, ++index)
      if (0 != propertiesE[index])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldsVertex[i]/propertiesE[index],
				     tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[index], fieldsVertex[i],
				     tolerance);
  } // for

  // Test vertex array sizes.
  size_t size = data.numPropsVertex + data.numVarsVertex;
  CPPUNIT_ASSERT_EQUAL(size, friction._propsStateVarsVertex.size());
} // testInitialize

// ----------------------------------------------------------------------
// Test getField().
void
pylith::friction::TestFrictionModel::testGetField(void)
{ // testGetField
  const double fieldE[] = { 0.6, 0.4 };
  const int fiberDim = 1;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const topology::Field<topology::SubMesh>& field = 
    friction.getField("friction_coefficient");

  const ALE::Obj<SieveSubMesh>& sieveMesh = field.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const double tolerance = 1.0e-06;

  int index = 0;
  double_array fieldVertex(fiberDim);
  const ALE::Obj<RealSection>& fieldSection = field.section();
  CPPUNIT_ASSERT(!fieldSection.isNull());
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldSection->getFiberDimension(*v_iter));
    fieldSection->restrictPoint(*v_iter, &fieldVertex[0], fieldVertex.size());
    for (int i=0; i < fiberDim; ++i, ++index)
      if (0 != fieldE[index])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldVertex[i]/fieldE[index], tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(fieldE[index], fieldVertex[i], tolerance);
  } // for
} // testGetField

// ----------------------------------------------------------------------
// Test retrievePropsStateVars().
void
pylith::friction::TestFrictionModel::testRetrievePropsStateVars(void)
{ // testRetrievePropsStateVars
  const double propertiesE[] = { 0.4, 1000000 };
  const size_t numProperties = 2;
  const double* stateVarsE = 0;
  const size_t numStateVars = 0;
  const int vertex = 2;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  friction.retrievePropsStateVars(vertex);

  const double tolerance = 1.0e-06;

  const double_array& fieldsVertex = friction._propsStateVarsVertex;
  CPPUNIT_ASSERT_EQUAL(numProperties + numStateVars, fieldsVertex.size());

  // Check properties array.
  int index = 0;
  CPPUNIT_ASSERT(0 != propertiesE);
  for (size_t i=0; i < numProperties; ++i)
    if (0.0 != propertiesE[i])
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldsVertex[index++]/propertiesE[i],
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], fieldsVertex[index++],
				   tolerance);

  // Check state variables array.
  CPPUNIT_ASSERT( (0 < numStateVars && 0 != stateVarsE) ||
		  (0 == numStateVars && 0 == stateVarsE) );
  for (size_t i=0; i < numStateVars; ++i)
    if (0.0 != stateVarsE[i])
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldsVertex[index++]/stateVarsE[i],
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], fieldsVertex[index++],
				   tolerance);
} // testRetrievePropsStateVars

// ----------------------------------------------------------------------
// Test calcFriction()
void
pylith::friction::TestFrictionModel::testCalcFriction(void)
{ // testCalcFriction
  const double slip = 1.2;
  const double slipRate = -2.3;
  const double normalTraction = -2.4;
  const double frictionCoef = 0.45;
  const double cohesion = 1000000;
  const double frictionE = -normalTraction*frictionCoef + cohesion;
  const int vertex = 2;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  friction.timeStep(data.dt);
  friction.retrievePropsStateVars(vertex);
  const double frictionV = friction.calcFriction(slip, slipRate, normalTraction);

  const double tolerance = 1.0e-6;
  if (0.0 != frictionE)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, frictionV/frictionE, tolerance);
  else
    CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionE, frictionV, tolerance);
} // testCalcFriction
    
// ----------------------------------------------------------------------
// Test updateStateVars()
void
pylith::friction::TestFrictionModel::testUpdateStateVars(void)
{ // testUpdateStateVars
  { // Test with friction model without state variables
    topology::Mesh mesh;
    faults::FaultCohesiveDyn fault;
    StaticFriction friction;
    StaticFrictionData data;
    _initialize(&mesh, &fault, &friction, &data);
    
    const double slip = 1.2;
    const double slipRate = -2.3;
    const double normalTraction = -2.4;
    const double cohesion = 1000000;
    const int vertex = 2;
    
    friction.timeStep(data.dt);
    friction.retrievePropsStateVars(vertex);
    friction.updateStateVars(slip, slipRate, normalTraction, vertex);
    
    // no outcome to test
  } // Test with friction model without state variables

  { // Test with friction model with state variables (slip weakening)
    // Initialize uses static friction, so we change to slip weakening.
    topology::Mesh mesh;
    faults::FaultCohesiveDyn fault;
    StaticFriction frictionDummy;
    StaticFrictionData data;
    _initialize(&mesh, &fault, &frictionDummy, &data);
    
    SlipWeakening friction;
    spatialdata::spatialdb::SimpleDB db;
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename("data/friction_slipweakening.spatialdb");
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
    friction.dbProperties(&db);
    fault.frictionModel(&friction);
    
    const double upDir[] = { 0.0, 0.0, 1.0 };
    fault.initialize(mesh, upDir);
    const int vertex = 2;

    const double slip = 0.25;
    const double slipRate = 0.64;
    const double normalTraction = -2.3;
    const double cohesion = 1000000;
    const double dt = 0.01;

    const int numStateVars = 2;
    const double stateVars[2] = { 0.5, 0.1 };
    const double stateVarsUpdatedE[2] = { 0.65, 0.25 };
    

    // Set state variables to given values
    friction._propsStateVarsVertex = 0.0;
    CPPUNIT_ASSERT_EQUAL(numStateVars, friction._varsFiberDim);
    for (size_t i=0; i < numStateVars; ++i)
      friction._propsStateVarsVertex[friction._propsFiberDim+i] = stateVars[i];

    friction.timeStep(dt);
    friction.updateStateVars(slip, slipRate, normalTraction, vertex);
    
    const double tolerance = 1.0e-06;
    CPPUNIT_ASSERT(0 != friction._fieldsPropsStateVars);
    const ALE::Obj<SubRealUniformSection>& fieldsSection = 
      friction._fieldsPropsStateVars->section();
    CPPUNIT_ASSERT_EQUAL(friction._fieldsPropsStateVars->fiberDim(),
			 fieldsSection->getFiberDimension(vertex));
    const double* fieldsVertex = fieldsSection->restrictPoint(vertex);
    CPPUNIT_ASSERT(fieldsVertex);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsUpdatedE[0],
				 fieldsVertex[friction._propsFiberDim+0], 
				 tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsUpdatedE[1],
				 fieldsVertex[friction._propsFiberDim+1],
				 tolerance);
  } // Test with friction model with state variables (slip weakening)

} // testUpdateStateVars

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestFrictionModel::setUp(void)
{ // setUp
  _friction = 0;
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::friction::TestFrictionModel::tearDown(void)
{ // tearDown
  delete _friction; _friction = 0;
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::friction::TestFrictionModel::testDBToProperties(void)
{ // testDBToProperties
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  // Check to make sure names of Metadata values match names of test
  // data values (consistency check).
  const int numDBProperties = _data->numDBProperties;
  char** dbPropertyLabelsE = _data->dbPropertyValues;
  CPPUNIT_ASSERT_EQUAL(numDBProperties, _friction->_metadata.numDBProperties());
  const char* const* dbPropertyLabels = _friction->_metadata.dbProperties();
  for (int i=0; i < numDBProperties; ++i) 
    CPPUNIT_ASSERT_EQUAL(std::string(dbPropertyLabelsE[i]),
			 std::string(dbPropertyLabels[i]));

  // Test _dbToProperties()
  const int numLocs = _data->numLocs;
  double_array dbValues(numDBProperties);

  const int propertiesSize = _data->numPropsVertex;
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBProperties; ++i)
      dbValues[i] = _data->dbProperties[iLoc*numDBProperties+i];

    _friction->_dbToProperties(&properties[0], dbValues);
    
    const double* const propertiesE = &_data->properties[iLoc*propertiesSize];
    const double tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     properties[i]/propertiesE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], properties[i],
				     tolerance);
    } // for
  } // for
} // testDBToProperties

// ----------------------------------------------------------------------
// Test _nondimProperties().
void
pylith::friction::TestFrictionModel::testNonDimProperties(void)
{ // testNonDimProperties
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsVertex;
  double_array propertiesNondim(propertiesSize);
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < propertiesSize; ++i)
      properties[i] = _data->properties[iLoc*propertiesSize+i];
    _friction->_nondimProperties(&properties[0], properties.size());
    
    const double* const propertiesNondimE =
      &_data->propertiesNondim[iLoc*propertiesSize];
    CPPUNIT_ASSERT(0 != propertiesNondimE);

    const double tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesNondimE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     properties[i]/propertiesNondimE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesNondimE[i], properties[i],
				     tolerance);
    } // for
  } // for
} // testNonDimProperties

// ----------------------------------------------------------------------
// Test _dimProperties().
void
pylith::friction::TestFrictionModel::testDimProperties(void)
{ // testDimProperties
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsVertex;
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < propertiesSize; ++i)
      properties[i] = _data->propertiesNondim[iLoc*propertiesSize+i];
    _friction->_dimProperties(&properties[0], properties.size());
    
    const double* const propertiesE =
      &_data->properties[iLoc*propertiesSize];
    CPPUNIT_ASSERT(0 != propertiesE);

    const double tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     properties[i]/propertiesE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], properties[i],
				     tolerance);
    } // for
  } // for
} // testDimProperties

// ----------------------------------------------------------------------
// Test _dbToStateVars().
void
pylith::friction::TestFrictionModel::testDBToStateVars(void)
{ // testDBToStateVars
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  // Check to make sure names of Metadata values match names of test
  // data values (consistency check).
  const int numDBStateVars = _data->numDBStateVars;
  char** dbStateVarsLabelsE = _data->dbStateVarValues;
  CPPUNIT_ASSERT_EQUAL(numDBStateVars, _friction->_metadata.numDBStateVars());
  const char* const* dbStateVarsLabels = _friction->_metadata.dbStateVars();
  for (int i=0; i < numDBStateVars; ++i) 
    CPPUNIT_ASSERT_EQUAL(std::string(dbStateVarsLabelsE[i]),
			 std::string(dbStateVarsLabels[i]));

  // Test _dbToStateVars()
  const int numLocs = _data->numLocs;
  double_array dbValues(numDBStateVars);

  const int stateVarsSize = _data->numVarsVertex;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBStateVars; ++i)
      dbValues[i] = _data->dbStateVars[iLoc*numDBStateVars+i];

    _friction->_dbToStateVars(&stateVars[0], dbValues);
    
    const double* const stateVarsE = 
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && 0 != stateVarsE) ||
		    (0 == stateVarsSize && 0 == stateVarsE) );
    const double tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     stateVars[i]/stateVarsE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i],
				     tolerance);
    } // for
  } // for
} // testDBToStateVars

// ----------------------------------------------------------------------
// Test _nondimStateVars().
void
pylith::friction::TestFrictionModel::testNonDimStateVars(void)
{ // testNonDimStateVars
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsVertex;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < stateVarsSize; ++i)
      stateVars[i] = _data->dbStateVars[iLoc*stateVarsSize+i];
    _friction->_nondimStateVars(&stateVars[0], stateVars.size());
    
    const double* const stateVarsNondimE =
      (stateVarsSize > 0) ? &_data->stateVarsNondim[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && 0 != stateVarsNondimE) ||
		    (0 == stateVarsSize && 0 == stateVarsNondimE) );

    const double tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsNondimE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     stateVars[i]/stateVarsNondimE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsNondimE[i], stateVars[i],
				     tolerance);
    } // for
  } // for
} // testNonDimStateVars

// ----------------------------------------------------------------------
// Test _dimStateVars().
void
pylith::friction::TestFrictionModel::testDimStateVars(void)
{ // testDimStateVars
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsVertex;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < stateVarsSize; ++i)
      stateVars[i] = _data->stateVarsNondim[iLoc*stateVarsSize+i];
    _friction->_dimStateVars(&stateVars[0], stateVars.size());
    
    const double* const stateVarsE =
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && 0 != stateVarsE) ||
		    (0 == stateVarsSize && 0 == stateVarsE) );

    const double tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     stateVars[i]/stateVarsE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i],
				     tolerance);
    } // for
  } // for
} // testDimStateVars

// ----------------------------------------------------------------------
// Test _calcFriction()
void
pylith::friction::TestFrictionModel::test_calcFriction(void)
{ // _testCalcFriction
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);

  const int numLocs = _data->numLocs;
  const int numPropsVertex = _data->numPropsVertex;
  const int numVarsVertex = _data->numVarsVertex;
  
  double_array properties(numPropsVertex);
  double_array stateVars(numVarsVertex);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numPropsVertex; ++i)
      properties[i] = _data->properties[iLoc*numPropsVertex+i];
    for (int i=0; i < numVarsVertex; ++i)
      stateVars[i] = _data->stateVars[iLoc*numVarsVertex+i];
    const double slip = _data->slip[iLoc];
    const double slipRate = _data->slipRate[iLoc];
    const double normalTraction = _data->normalTraction[iLoc];

    _friction->timeStep(_data->dt);
    const double friction = _friction->_calcFriction(
					slip, slipRate, normalTraction,
					&properties[0], properties.size(),
					&stateVars[0], stateVars.size());
    
    const double frictionE = _data->friction[iLoc];
    
    const double tolerance = 1.0e-06;

    if (0.0 != frictionE)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, friction/frictionE, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionE, friction, tolerance);
  } // for
} // _testCalcFriction

// ----------------------------------------------------------------------
// Test _updateStateVars()
void
pylith::friction::TestFrictionModel::test_updateStateVars(void)
{ // test_updateStateVars
  CPPUNIT_ASSERT(0 != _friction);
  CPPUNIT_ASSERT(0 != _data);

  const bool computeStateVars = true;

  const int numLocs = _data->numLocs;
  const int numPropsVertex = _data->numPropsVertex;
  const int numVarsVertex = _data->numVarsVertex;
  
  double_array properties(numPropsVertex);
  double_array stateVars(numVarsVertex);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    const double slip = _data->slip[iLoc];
    const double slipRate = _data->slipRate[iLoc];
    const double normalTraction = _data->normalTraction[iLoc];
    for (int i=0; i < numPropsVertex; ++i)
      properties[i] = _data->properties[iLoc*numPropsVertex+i];
    for (int i=0; i < numVarsVertex; ++i)
      stateVars[i] = _data->stateVars[iLoc*numVarsVertex+i];

    _friction->timeStep(_data->dt);
    _friction->_updateStateVars(slip, slipRate, normalTraction,
				&stateVars[0], stateVars.size(),
				&properties[0], properties.size());
    
    const double* stateVarsE = 
      (numVarsVertex > 0) ? &_data->stateVarsUpdated[iLoc*numVarsVertex] : 0;
    CPPUNIT_ASSERT( (0 < numVarsVertex && 0 != stateVarsE) ||
		    (0 == numVarsVertex && 0 == stateVarsE) );

    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsVertex; ++i) {
#if 1 // DEBUGGING
      std::cout << "valE: " << stateVarsE[i] 
		<< ", val: " << stateVars[i]
		<< std::endl;
#endif
      if (0.0 != stateVarsE[i])
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i],
				     tolerance);
    } // for
  } // for
} // test_updateStateVars

// ----------------------------------------------------------------------
// Setup nondimensionalization.
void
pylith::friction::TestFrictionModel::setupNormalizer(void)
{ // setupNormalizer
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.timeScale(_data->timeScale);
  normalizer.densityScale(_data->densityScale);
  _friction->normalizer(normalizer);
} // setupNormalizer

// ----------------------------------------------------------------------
// Setup mesh and material.
void
pylith::friction::TestFrictionModel::_initialize(
					  topology::Mesh* mesh,
					  faults::FaultCohesiveDyn* fault,
					  StaticFriction* friction,
					  const StaticFrictionData* data)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != friction);
  CPPUNIT_ASSERT(0 != data);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  // Setup quadrature
  feassemble::Quadrature<topology::SubMesh> quadrature;
  feassemble::GeometryLine2D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const double basis[4] = { 1.0, 0.0, 0.0, 1.0 };
  const double basisDeriv[4] = { -0.5, 0.5, -0.5, 0.5 };
  const double quadPtsRef[2] = { -1.0, 1.0 };
  const double quadWts[2] = { 1.0, 1.0  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  const bool flipFault = false;
  const char* label = "fault";
  int firstFaultVertex = 0;
  int firstLagrangeVertex = mesh->sieveMesh()->getIntSection(label)->size();
  int firstFaultCell = mesh->sieveMesh()->getIntSection(label)->size();
  if (fault->useLagrangeConstraints())
    firstFaultCell += mesh->sieveMesh()->getIntSection(label)->size();
  fault->id(100);
  fault->label(label);
  fault->quadrature(&quadrature);
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex,
      &firstFaultCell, flipFault);

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename("data/friction_static.spatialdb");
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  friction->dbProperties(&db);
  friction->label("my_friction");
  friction->normalizer(normalizer);
  fault->frictionModel(friction);
  
  const double upDir[3] = { 0.0, 0.0, 1.0 };
  fault->initialize(*mesh, upDir);
  fault->verifyConfiguration(*mesh);
} // _initialize


// End of file 
