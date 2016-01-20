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

#include "TestFrictionModel.hh" // Implementation of class methods

#include "data/StaticFrictionData.hh" // USES StaticFrictionData
#include "data/SlipWeakeningData.hh" // USES SlipWeakeningData

#include "pylith/friction/StaticFriction.hh" // USES StaticFriction
#include "pylith/friction/SlipWeakening.hh" // USES SlipWeakening

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine2D.hh" // USES GeometryLine2D

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestFrictionModel );

// ----------------------------------------------------------------------
// Test label()
void
pylith::friction::TestFrictionModel::testLabel(void)
{ // testLabel
  PYLITH_METHOD_BEGIN;

  const std::string& label = "the database";
  StaticFriction friction;
  friction.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction.label()));

  PYLITH_METHOD_END;
} // testLabel
    
// ----------------------------------------------------------------------
// Test timestep()
void
pylith::friction::TestFrictionModel::testTimeStep(void) 
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  const PylithScalar dt = 2.0;
  StaticFriction friction;
  friction.timeStep(dt);
  
  CPPUNIT_ASSERT_EQUAL(dt, friction.timeStep());

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test dbProperties()
void
pylith::friction::TestFrictionModel::testDBProperties(void)
{ // testDBProperties
  PYLITH_METHOD_BEGIN;

  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  StaticFriction friction;
  friction.dbProperties(&db);
  
  CPPUNIT_ASSERT(friction._dbProperties);
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction._dbProperties->label()));

  PYLITH_METHOD_END;
} // testDBProperties

// ----------------------------------------------------------------------
// Test dbStateVars()
void
pylith::friction::TestFrictionModel::testDBStateVars(void)
{ // testDBStateVars
  PYLITH_METHOD_BEGIN;

  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  StaticFriction friction;
  friction.dbInitialState(&db);
  
  CPPUNIT_ASSERT(friction._dbInitialState);
  CPPUNIT_ASSERT_EQUAL(label, std::string(friction._dbInitialState->label()));

  PYLITH_METHOD_END;
} // testDBStateVars

// ----------------------------------------------------------------------
// Test normalizer()
void
pylith::friction::TestFrictionModel::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;

  spatialdata::units::Nondimensional normalizer;
  const double lengthScale = 2.0;
  normalizer.lengthScale(lengthScale);

  StaticFriction friction;
  friction.normalizer(normalizer);
  
  CPPUNIT_ASSERT(friction._normalizer);
  CPPUNIT_ASSERT_EQUAL(lengthScale, friction._normalizer->lengthScale());

  PYLITH_METHOD_END;
} // testNormalizer

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::friction::TestFrictionModel::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);
  CPPUNIT_ASSERT(friction._fieldsPropsStateVars);

  const PylithScalar propertiesE[2*2] = {
    0.6, 1000000/data.pressureScale,
    0.4, 1000000/data.pressureScale,
  };

  PetscDM faultDMMesh = fault.faultMesh().dmMesh();CPPUNIT_ASSERT(faultDMMesh);
  topology::Stratum depthStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const PylithScalar tolerance = 1.0e-06;

  // Test fieldsPropsStateVars with mesh
  const materials::Metadata& metadata = friction.getMetadata();
  const int numProperties = metadata.numProperties();
  int index = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    for (int i = 0; i < numProperties; ++i, ++index) {
      const materials::Metadata::ParamDescription& property = metadata.getProperty(i);
      topology::Field& prop = friction._fieldsPropsStateVars->get(property.name.c_str());
      topology::VecVisitorMesh propVisitor(prop);
      const PetscScalar* propArray = propVisitor.localArray();CPPUNIT_ASSERT(propArray);

      const PetscInt off = propVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(1, propVisitor.sectionDof(v));
      if (0.0 != propertiesE[index]) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propArray[off]/propertiesE[index], tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[index], propArray[off], tolerance);
      } // if/else
    } // for
  } // for

  // Test vertex array sizes.
  size_t size = data.numPropsVertex + data.numVarsVertex;
  CPPUNIT_ASSERT_EQUAL(size, friction._propsStateVarsVertex.size());

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test getField().
void
pylith::friction::TestFrictionModel::testGetField(void)
{ // testGetField
  PYLITH_METHOD_BEGIN;

  const PylithScalar fieldE[] = { 0.6, 0.4 };
  const int fiberDim = 1;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const topology::Field& frictionField = friction.getField("friction_coefficient");
  topology::VecVisitorMesh frictionVisitor(frictionField);
  PetscScalar *frictionArray = frictionVisitor.localArray();CPPUNIT_ASSERT(frictionArray);

  PetscDM dmMesh = frictionField.mesh().dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum depthStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const PylithScalar tolerance = 1.0e-06;
  int index = 0;

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = frictionVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, frictionVisitor.sectionDof(v));
    for(int i = 0; i < fiberDim; ++i, ++index)
      if (0.0 != fieldE[index])
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, frictionArray[off+i]/fieldE[index], tolerance);
      else
        CPPUNIT_ASSERT_DOUBLES_EQUAL(fieldE[index], frictionArray[off+i], tolerance);
  } // for

  PYLITH_METHOD_END;
} // testGetField

// ----------------------------------------------------------------------
// Test retrievePropsStateVars().
void
pylith::friction::TestFrictionModel::testRetrievePropsStateVars(void)
{ // testRetrievePropsStateVars
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const size_t numProperties = 2;
  const PylithScalar propertiesE[numProperties] = {
    0.6, 1000000/data.pressureScale,
  };
  const PylithScalar* stateVarsE = 0;
  const size_t numStateVars = 0;
  const int vertex = 2;

  friction.retrievePropsStateVars(vertex);

  const PylithScalar tolerance = 1.0e-06;

  const scalar_array& fieldsVertex = friction._propsStateVarsVertex;
  CPPUNIT_ASSERT_EQUAL(numProperties + numStateVars, fieldsVertex.size());

  // Check properties array.
  int index = 0;
  CPPUNIT_ASSERT(propertiesE);
  for (size_t i=0; i < numProperties; ++i) {
    std::cout << "propertiesE: " << propertiesE[i] << ", properties: " << fieldsVertex[index] << std::endl;
    if (0.0 != propertiesE[i])
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldsVertex[index++]/propertiesE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], fieldsVertex[index++], tolerance);
  } // for

  // Check state variables array.
  CPPUNIT_ASSERT( (0 < numStateVars && stateVarsE) ||
		  (0 == numStateVars && !stateVarsE) );
  for (size_t i=0; i < numStateVars; ++i)
    if (0.0 != stateVarsE[i])
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldsVertex[index++]/stateVarsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], fieldsVertex[index++], tolerance);

  PYLITH_METHOD_END;
} // testRetrievePropsStateVars

// ----------------------------------------------------------------------
// Test calcFriction()
void
pylith::friction::TestFrictionModel::testCalcFriction(void)
{ // testCalcFriction
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const PylithScalar t = 1.5;
  const PylithScalar slip = 1.2;
  const PylithScalar slipRate = -2.3;
  const PylithScalar normalTraction = -2.4e-3;
  const PylithScalar frictionCoef = 0.6;
  const PylithScalar cohesion = 1.0e+6/data.pressureScale;
  const PylithScalar frictionE = -normalTraction*frictionCoef + cohesion;
  const int vertex = 2;

  friction.timeStep(data.dt);
  friction.retrievePropsStateVars(vertex);
  const PylithScalar frictionV = friction.calcFriction(t, slip, slipRate, normalTraction);

  const PylithScalar tolerance = 1.0e-6;
  if (0.0 != frictionE) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, frictionV/frictionE, tolerance);
  } else {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionE, frictionV, tolerance);
  } // if/else

  PYLITH_METHOD_END;
} // testCalcFriction
    
// ----------------------------------------------------------------------
// Test calcFrictionDeriv()
void
pylith::friction::TestFrictionModel::testCalcFrictionDeriv(void)
{ // testCalcFrictionDeriv
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  faults::FaultCohesiveDyn fault;
  StaticFriction friction;
  StaticFrictionData data;
  _initialize(&mesh, &fault, &friction, &data);

  const PylithScalar t = 1.5;
  const PylithScalar slip = 1.2;
  const PylithScalar slipRate = -2.3;
  const PylithScalar normalTraction = -2.4e-3;
  const PylithScalar frictionCoef = 0.6;
  const PylithScalar cohesion = 1.0e+6/data.pressureScale;
  const PylithScalar frictionDerivE = 0.0;
  const int vertex = 2;

  friction.timeStep(data.dt);
  friction.retrievePropsStateVars(vertex);
  const PylithScalar frictionDeriv = friction.calcFrictionDeriv(t, slip, slipRate, normalTraction);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionDerivE, frictionDeriv, tolerance);

  PYLITH_METHOD_END;
} // testCalcFrictionDeriv
    
// ----------------------------------------------------------------------
// Test updateStateVars()
void
pylith::friction::TestFrictionModel::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  { // Test with friction model without state variables
    topology::Mesh mesh;
    faults::FaultCohesiveDyn fault;
    StaticFriction friction;
    StaticFrictionData data;
    _initialize(&mesh, &fault, &friction, &data);
    
    const PylithScalar t = 1.5;
    const PylithScalar slip = 1.2;
    const PylithScalar slipRate = -2.3;
    const PylithScalar normalTraction = -2.4;
    const PylithScalar cohesion = 1000000;
    const int vertex = 2;
    
    friction.timeStep(data.dt);
    friction.retrievePropsStateVars(vertex);
    friction.updateStateVars(t, slip, slipRate, normalTraction, vertex);
    
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
    
    const PylithScalar upDir[] = { 0.0, 0.0, 1.0 };
    fault.initialize(mesh, upDir);
    const int vertex = 2;

    const PylithScalar t = 1.5;
    const PylithScalar slip = 0.25;
    const PylithScalar slipRate = 0.64;
    const PylithScalar normalTraction = -2.3;
    const PylithScalar cohesion = 1000000;
    const PylithScalar dt = 0.01;

    const PylithScalar stateVars[2] = { 0.5, 0.1 };
    const PylithScalar stateVarsUpdatedE[2] = { 0.65, 0.25 };
    
    const materials::Metadata& metadata = friction.getMetadata();
    const int numStateVars = metadata.numStateVars();
    CPPUNIT_ASSERT(2 == numStateVars);

    // Set state variables to given values
    friction._propsStateVarsVertex = 0.0;
    CPPUNIT_ASSERT_EQUAL(numStateVars, friction._varsFiberDim);
    for (size_t i=0; i < numStateVars; ++i)
      friction._propsStateVarsVertex[friction._propsFiberDim+i] = stateVars[i];

    friction.timeStep(dt);
    friction.updateStateVars(t, slip, slipRate, normalTraction, vertex);
    
    const PylithScalar tolerance = 1.0e-06;
    CPPUNIT_ASSERT(friction._fieldsPropsStateVars);
    for(PetscInt i = 0; i < numStateVars; ++i) {
      const materials::Metadata::ParamDescription& stateVar = metadata.getStateVar(i);
      topology::Field& stateVarField = friction._fieldsPropsStateVars->get(stateVar.name.c_str());
      topology::VecVisitorMesh stateVarVisitor(stateVarField);
      PetscScalar *fieldsArray = stateVarVisitor.localArray();CPPUNIT_ASSERT(fieldsArray);

      const PetscInt off = stateVarVisitor.sectionOffset(vertex);
      CPPUNIT_ASSERT_EQUAL(1, stateVarVisitor.sectionDof(vertex));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsUpdatedE[i], fieldsArray[off], tolerance);
    } // for
  } // Test with friction model with state variables (slip weakening)

  PYLITH_METHOD_END;
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestFrictionModel::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _friction = 0;
  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::friction::TestFrictionModel::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _friction; _friction = 0;
  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::friction::TestFrictionModel::testDBToProperties(void)
{ // testDBToProperties
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
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
  scalar_array dbValues(numDBProperties);

  const int propertiesSize = _data->numPropsVertex;
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBProperties; ++i)
      dbValues[i] = _data->dbProperties[iLoc*numDBProperties+i];

    _friction->_dbToProperties(&properties[0], dbValues);
    
    const PylithScalar* const propertiesE = &_data->properties[iLoc*propertiesSize];
    const PylithScalar tolerance = 1.0e-06;
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

  PYLITH_METHOD_END;
} // testDBToProperties

// ----------------------------------------------------------------------
// Test _nondimProperties().
void
pylith::friction::TestFrictionModel::testNonDimProperties(void)
{ // testNonDimProperties
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsVertex;
  scalar_array propertiesNondim(propertiesSize);
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < propertiesSize; ++i)
      properties[i] = _data->properties[iLoc*propertiesSize+i];
    _friction->_nondimProperties(&properties[0], properties.size());
    
    const PylithScalar* const propertiesNondimE =
      &_data->propertiesNondim[iLoc*propertiesSize];
    CPPUNIT_ASSERT(propertiesNondimE);

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesNondimE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, properties[i]/propertiesNondimE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesNondimE[i], properties[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testNonDimProperties

// ----------------------------------------------------------------------
// Test _dimProperties().
void
pylith::friction::TestFrictionModel::testDimProperties(void)
{ // testDimProperties
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsVertex;
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < propertiesSize; ++i)
      properties[i] = _data->propertiesNondim[iLoc*propertiesSize+i];
    _friction->_dimProperties(&properties[0], properties.size());
    
    const PylithScalar* const propertiesE = &_data->properties[iLoc*propertiesSize];
    CPPUNIT_ASSERT(propertiesE);

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, properties[i]/propertiesE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], properties[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testDimProperties

// ----------------------------------------------------------------------
// Test _dbToStateVars().
void
pylith::friction::TestFrictionModel::testDBToStateVars(void)
{ // testDBToStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
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
  scalar_array dbValues(numDBStateVars);

  const int stateVarsSize = _data->numVarsVertex;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBStateVars; ++i)
      dbValues[i] = _data->dbStateVars[iLoc*numDBStateVars+i];

    _friction->_dbToStateVars(&stateVars[0], dbValues);
    
    const PylithScalar* const stateVarsE = 
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsE) ||
		    (0 == stateVarsSize && !stateVarsE) );
    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testDBToStateVars

// ----------------------------------------------------------------------
// Test _nondimStateVars().
void
pylith::friction::TestFrictionModel::testNonDimStateVars(void)
{ // testNonDimStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsVertex;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < stateVarsSize; ++i)
      stateVars[i] = _data->dbStateVars[iLoc*stateVarsSize+i];
    _friction->_nondimStateVars(&stateVars[0], stateVars.size());
    
    const PylithScalar* const stateVarsNondimE =
      (stateVarsSize > 0) ? &_data->stateVarsNondim[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsNondimE) ||
		    (0 == stateVarsSize && !stateVarsNondimE) );

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsNondimE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsNondimE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsNondimE[i], stateVars[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testNonDimStateVars

// ----------------------------------------------------------------------
// Test _dimStateVars().
void
pylith::friction::TestFrictionModel::testDimStateVars(void)
{ // testDimStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsVertex;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < stateVarsSize; ++i)
      stateVars[i] = _data->stateVarsNondim[iLoc*stateVarsSize+i];
    _friction->_dimStateVars(&stateVars[0], stateVars.size());
    
    const PylithScalar* const stateVarsE =
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsE) ||
		    (0 == stateVarsSize && !stateVarsE) );

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testDimStateVars

// ----------------------------------------------------------------------
// Test _calcFriction()
void
pylith::friction::TestFrictionModel::test_calcFriction(void)
{ // _testCalcFriction
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);

  const int numLocs = _data->numLocs;
  const int numPropsVertex = _data->numPropsVertex;
  const int numVarsVertex = _data->numVarsVertex;
  
  scalar_array properties(numPropsVertex);
  scalar_array stateVars(numVarsVertex);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numPropsVertex; ++i)
      properties[i] = _data->properties[iLoc*numPropsVertex+i];
    for (int i=0; i < numVarsVertex; ++i)
      stateVars[i] = _data->stateVars[iLoc*numVarsVertex+i];
    const PylithScalar t = 1.5;
    const PylithScalar slip = _data->slip[iLoc];
    const PylithScalar slipRate = _data->slipRate[iLoc];
    const PylithScalar normalTraction = _data->normalTraction[iLoc];

    _friction->timeStep(_data->dt);
    const PylithScalar friction = 
      _friction->_calcFriction(t, slip, slipRate, normalTraction,
			       &properties[0], properties.size(),
			       &stateVars[0], stateVars.size());
    
    const PylithScalar frictionE = _data->friction[iLoc];
    
    const PylithScalar tolerance = 1.0e-06;

    if (0.0 != frictionE)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, friction/frictionE, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionE, friction, tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testCalcFriction

// ----------------------------------------------------------------------
// Test _calcFriction()
void
pylith::friction::TestFrictionModel::test_calcFrictionDeriv(void)
{ // _testCalcFrictionDeriv
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);

  const int numLocs = _data->numLocs;
  const int numPropsVertex = _data->numPropsVertex;
  const int numVarsVertex = _data->numVarsVertex;
  
  scalar_array properties(numPropsVertex);
  scalar_array stateVars(numVarsVertex);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numPropsVertex; ++i)
      properties[i] = _data->properties[iLoc*numPropsVertex+i];
    for (int i=0; i < numVarsVertex; ++i)
      stateVars[i] = _data->stateVars[iLoc*numVarsVertex+i];
    const PylithScalar t = 1.5;
    const PylithScalar slip = _data->slip[iLoc];
    const PylithScalar slipRate = _data->slipRate[iLoc];
    const PylithScalar normalTraction = _data->normalTraction[iLoc];

    _friction->timeStep(_data->dt);
    const PylithScalar frictionDeriv = 
      _friction->_calcFrictionDeriv(t, slip, slipRate, normalTraction,
				    &properties[0], properties.size(),
				    &stateVars[0], stateVars.size());
    
    const PylithScalar frictionDerivE = _data->frictionDeriv[iLoc];
    
    const PylithScalar tolerance = 1.0e-06;
    if (0.0 != frictionDerivE)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, frictionDeriv/frictionDerivE, tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(frictionDerivE, frictionDeriv, tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testCalcFrictionDeriv

// ----------------------------------------------------------------------
// Test _updateStateVars()
void
pylith::friction::TestFrictionModel::test_updateStateVars(void)
{ // test_updateStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_friction);
  CPPUNIT_ASSERT(_data);

  const bool computeStateVars = true;

  const int numLocs = _data->numLocs;
  const int numPropsVertex = _data->numPropsVertex;
  const int numVarsVertex = _data->numVarsVertex;
  
  scalar_array properties(numPropsVertex);
  scalar_array stateVars(numVarsVertex);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    const PylithScalar t = 1.5;
    const PylithScalar slip = _data->slip[iLoc];
    const PylithScalar slipRate = _data->slipRate[iLoc];
    const PylithScalar normalTraction = _data->normalTraction[iLoc];
    for (int i=0; i < numPropsVertex; ++i)
      properties[i] = _data->properties[iLoc*numPropsVertex+i];
    for (int i=0; i < numVarsVertex; ++i)
      stateVars[i] = _data->stateVars[iLoc*numVarsVertex+i];

    _friction->timeStep(_data->dt);
    _friction->_updateStateVars(t, slip, slipRate, normalTraction,
				&stateVars[0], stateVars.size(),
				&properties[0], properties.size());
    
    const PylithScalar* stateVarsE = 
      (numVarsVertex > 0) ? &_data->stateVarsUpdated[iLoc*numVarsVertex] : 0;
    CPPUNIT_ASSERT( (0 < numVarsVertex && stateVarsE) ||
		    (0 == numVarsVertex && !stateVarsE) );

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < numVarsVertex; ++i) {
#if 0 // DEBUGGING
      std::cout << "valE: " << stateVarsE[i] 
		<< ", val: " << stateVars[i]
		<< std::endl;
#endif
      if (0.0 != stateVarsE[i])
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // test_updateStateVars

// ----------------------------------------------------------------------
// Setup nondimensionalization.
void
pylith::friction::TestFrictionModel::setupNormalizer(void)
{ // setupNormalizer
  PYLITH_METHOD_BEGIN;

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.timeScale(_data->timeScale);
  normalizer.densityScale(_data->densityScale);
  _friction->normalizer(normalizer);

  PYLITH_METHOD_END;
} // setupNormalizer

// ----------------------------------------------------------------------
// Setup mesh and material.
void
pylith::friction::TestFrictionModel::_initialize(topology::Mesh* mesh,
						 faults::FaultCohesiveDyn* fault,
						 StaticFriction* friction,
						 const StaticFrictionData* data)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(friction);
  CPPUNIT_ASSERT(data);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  // Setup coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  
  // Setup scales
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(data->lengthScale);
  normalizer.pressureScale(data->pressureScale);
  normalizer.timeScale(data->timeScale);
  normalizer.densityScale(data->densityScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Setup quadrature
  feassemble::Quadrature quadrature;
  feassemble::GeometryLine2D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const PylithScalar basis[4] = { 1.0, 0.0, 0.0, 1.0 };
  const PylithScalar basisDeriv[4] = { -0.5, 0.5, -0.5, 0.5 };
  const PylithScalar quadPtsRef[2] = { -1.0, 1.0 };
  const PylithScalar quadWts[2] = { 1.0, 1.0  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  const char* label = "fault";

  PetscInt labelSize;
  PetscErrorCode err = DMGetStratumSize(mesh->dmMesh(), label, 1, &labelSize);PYLITH_CHECK_ERROR(err);

  PetscInt firstFaultVertex    = 0;
  PetscInt firstLagrangeVertex = labelSize;
  PetscInt firstFaultCell      = labelSize;
  if (fault->useLagrangeConstraints())
    firstFaultCell += labelSize;
  fault->id(100);
  fault->label(label);
  fault->quadrature(&quadrature);
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  fault->normalizer(normalizer);

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename("data/friction_static.spatialdb");
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  friction->dbProperties(&db);
  friction->label("my_friction");
  friction->normalizer(normalizer);
  fault->frictionModel(friction);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  fault->initialize(*mesh, upDir);
  fault->verifyConfiguration(*mesh);

  PYLITH_METHOD_END;
} // _initialize


// End of file 
