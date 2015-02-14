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

#include "TestMaterial.hh" // Implementation of class methods

#include "data/MaterialData.hh" // USES MaterialData

#include "pylith/utils/array.hh" // USES scalar_array

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcmp()
#include <cassert> // USES assert()
#include <cmath> // USES sqrt()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaterial );

// ----------------------------------------------------------------------
// Test id()
void
pylith::materials::TestMaterial::testID(void)
{ // testID
  PYLITH_METHOD_BEGIN;
 
  const int id = 346;
  ElasticPlaneStrain material;
  material.id(id);
  
  CPPUNIT_ASSERT_EQUAL(id,  material.id());

  PYLITH_METHOD_END;
} // testID

// ----------------------------------------------------------------------
// Test label()
void
pylith::materials::TestMaterial::testLabel(void)
{ // testLabel
  PYLITH_METHOD_BEGIN;
 
  const std::string& label = "the database";
  ElasticPlaneStrain material;
  material.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(material.label()));

  PYLITH_METHOD_END;
} // testLabel
    
// ----------------------------------------------------------------------
// Test timestep()
void
pylith::materials::TestMaterial::testTimeStep(void) 
{ // testTimeStep
  PYLITH_METHOD_BEGIN;
 
  const PylithScalar dt = 2.0;
  ElasticPlaneStrain material;
  material.timeStep(dt);
  
  CPPUNIT_ASSERT_EQUAL(dt, material.timeStep());

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test dbProperties()
void
pylith::materials::TestMaterial::testDBProperties(void)
{ // testDBProperties
  PYLITH_METHOD_BEGIN;
 
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticPlaneStrain material;
  material.dbProperties(&db);
  
  CPPUNIT_ASSERT(material._dbProperties);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbProperties->label()));

  PYLITH_METHOD_END;
} // testDBProperties

// ----------------------------------------------------------------------
// Test dbStateVars()
void
pylith::materials::TestMaterial::testDBStateVars(void)
{ // testDBStateVars
  PYLITH_METHOD_BEGIN;
 
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticPlaneStrain material;
  material.dbInitialState(&db);
  
  CPPUNIT_ASSERT(material._dbInitialState);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialState->label()));

  PYLITH_METHOD_END;
} // testDBStateVars

// ----------------------------------------------------------------------
// Test normalizer()
void
pylith::materials::TestMaterial::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;
 
  spatialdata::units::Nondimensional normalizer;
  const double lengthScale = 2.0;
  normalizer.lengthScale(lengthScale);

  ElasticPlaneStrain material;
  material.normalizer(normalizer);
  
  CPPUNIT_ASSERT(material._normalizer);
  CPPUNIT_ASSERT_EQUAL(lengthScale, material._normalizer->lengthScale());

  PYLITH_METHOD_END;
} // testNormalizer

// ----------------------------------------------------------------------
// Test needNewJacobian()
void
pylith::materials::TestMaterial::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  PYLITH_METHOD_BEGIN;
 
  ElasticPlaneStrain material;

  bool flag = false;
  material._needNewJacobian = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.needNewJacobian());

  flag = true;
  material._needNewJacobian = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.needNewJacobian());

  PYLITH_METHOD_END;
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test isJacobianSymmetric()
void
pylith::materials::TestMaterial::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  PYLITH_METHOD_BEGIN;
 
  ElasticPlaneStrain material;

  CPPUNIT_ASSERT_EQUAL(true, material.isJacobianSymmetric());

  bool flag = false;
  material._isJacobianSymmetric = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.isJacobianSymmetric());

  PYLITH_METHOD_END;
} // testIsJacobianSymmetric

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestMaterial::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;
 
  // Setup mesh
  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  const PylithScalar lengthScale = 1.0e+3;
  const PylithScalar pressureScale = 2.25e+10;
  const PylithScalar timeScale = 2.0;
  const PylithScalar velocityScale = lengthScale / timeScale;
  const PylithScalar densityScale = pressureScale / (velocityScale*velocityScale);
  normalizer.lengthScale(lengthScale);
  normalizer.pressureScale(pressureScale);
  normalizer.timeScale(timeScale);
  normalizer.densityScale(densityScale);
  topology::MeshOps::nondimensionalize(&mesh, normalizer);

  // Setup quadrature
  feassemble::Quadrature quadrature;
  feassemble::GeometryTri2D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 2;
  const int numCorners = 3;
  const int numQuadPts = 1;
  const int spaceDim = 2;
  const PylithScalar basis[] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
  const PylithScalar basisDeriv[] = { 
    -0.5, 0.5,
    -0.5, 0.0,
     0.0, 0.5,
  };
  const PylithScalar quadPtsRef[] = { -1.0/3.0, -1.0/3.0 };
  const PylithScalar quadWts[] = { 2.0  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);


  // Get cells associated with material
  const PetscInt  materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();

  // Compute geometry for cells
  quadrature.initializeGeometry();

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename("data/matinitialize.spatialdb");
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  ElasticPlaneStrain material;
  material.dbProperties(&db);
  material.id(materialId);
  material.label("my_material");
  material.normalizer(normalizer);
  material.initialize(mesh, &quadrature);

  const PylithScalar densityA = 2500.0;
  const PylithScalar vsA = 3000.0;
  const PylithScalar vpA = vsA*sqrt(3.0);
  const PylithScalar muA = vsA*vsA*densityA;
  const PylithScalar lambdaA = vpA*vpA*densityA - 2.0*muA;
  const PylithScalar densityB = 2000.0;
  const PylithScalar vsB = 1200.0;
  const PylithScalar vpB = vsB*sqrt(3.0);
  const PylithScalar muB = vsB*vsB*densityB;
  const PylithScalar lambdaB = vpB*vpB*densityB - 2.0*muB;
  const PylithScalar densityE[] = { densityA, densityB };
  const PylithScalar muE[] = { muA, muB };
  const PylithScalar lambdaE[] = { lambdaA, lambdaB };

  PetscInt cell = cells[0];
  const PylithScalar tolerance = 1.0e-06;

  CPPUNIT_ASSERT(material._properties);
  topology::VecVisitorMesh propertiesVisitor(*material._properties);
  const PetscScalar* propertiesArray = propertiesVisitor.localArray();
  const PetscInt off = propertiesVisitor.sectionOffset(cell);

  const int p_density = 0;
  const int p_mu = 1;
  const int p_lambda = 2;

  // density
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_density;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesArray[off+index]/densityE[i]*densityScale, tolerance);
  } // for
  
  // mu
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_mu;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesArray[off+index]/muE[i]*pressureScale, tolerance);
  } // for
  
  // lambda
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_lambda;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesArray[off+index]/lambdaE[i]*pressureScale, tolerance);
  } // for

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaterial::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;
 
  _material = 0;
  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::materials::TestMaterial::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;
 
  delete _material; _material = 0;
  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test dimension()
void
pylith::materials::TestMaterial::testDimension(void)
{ // testDimension
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT_EQUAL(_data->dimension, _material->dimension());

  PYLITH_METHOD_END;
} // testDimension

// ----------------------------------------------------------------------
// Test dimension()
void
pylith::materials::TestMaterial::testTensorSize(void)
{ // testTensorSize
  PYLITH_METHOD_BEGIN;
 
  int tensorSize = 0;
  const int dimension = _data->dimension;
  switch(dimension)
    { // switch
    case 1 :
      tensorSize = 1;
      break;
    case 2 :
      tensorSize = 3;
      break;
    case 3 :
      tensorSize = 6;
      break;
    default :
      assert(0);
    } // switch
  CPPUNIT_ASSERT(tensorSize > 0);

  CPPUNIT_ASSERT_EQUAL(tensorSize, _material->tensorSize());

  PYLITH_METHOD_END;
} // testTensorSize

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::materials::TestMaterial::testDBToProperties(void)
{ // testDBToProperties
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  // Check to make sure names of Metadata values match names of test
  // data values (consistency check).
  const int numDBProperties = _data->numDBProperties;
  char** dbPropertyLabelsE = _data->dbPropertyValues;
  CPPUNIT_ASSERT_EQUAL(numDBProperties, _material->_metadata.numDBProperties());
  const char* const* dbPropertyLabels = _material->_metadata.dbProperties();
  for (int i=0; i < numDBProperties; ++i) 
    CPPUNIT_ASSERT_EQUAL(std::string(dbPropertyLabelsE[i]),
			 std::string(dbPropertyLabels[i]));

  // Test _dbToProperties()
  const int numLocs = _data->numLocs;
  scalar_array dbValues(numDBProperties);

  const int propertiesSize = _data->numPropsQuadPt;
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBProperties; ++i)
      dbValues[i] = _data->dbProperties[iLoc*numDBProperties+i];

    _material->_dbToProperties(&properties[0], dbValues);
    
    const PylithScalar* const propertiesE = &_data->properties[iLoc*propertiesSize];
    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < propertiesSize; ++i) {
      if (fabs(propertiesE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, properties[i]/propertiesE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(propertiesE[i], properties[i], tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testDBToProperties

// ----------------------------------------------------------------------
// Test _nondimProperties().
void
pylith::materials::TestMaterial::testNonDimProperties(void)
{ // testNonDimProperties
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsQuadPt;
  scalar_array propertiesNondim(propertiesSize);
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&properties[0], &_data->properties[iLoc*propertiesSize],
	   propertiesSize*sizeof(PylithScalar));
    _material->_nondimProperties(&properties[0], properties.size());
    
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
pylith::materials::TestMaterial::testDimProperties(void)
{ // testDimProperties
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsQuadPt;
  scalar_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&properties[0], &_data->propertiesNondim[iLoc*propertiesSize], 
	   propertiesSize*sizeof(PylithScalar));
    _material->_dimProperties(&properties[0], properties.size());
    
    const PylithScalar* const propertiesE =
      &_data->properties[iLoc*propertiesSize];
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
pylith::materials::TestMaterial::testDBToStateVars(void)
{ // testDBToStateVars
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  // Check to make sure names of Metadata values match names of test
  // data values (consistency check).
  const int numDBStateVars = _data->numDBStateVars;
  char** dbStateVarsLabelsE = _data->dbStateVarValues;
  CPPUNIT_ASSERT_EQUAL(numDBStateVars, _material->_metadata.numDBStateVars());
  const char* const* dbStateVarsLabels = _material->_metadata.dbStateVars();
  for (int i=0; i < numDBStateVars; ++i) 
    CPPUNIT_ASSERT_EQUAL(std::string(dbStateVarsLabelsE[i]),
			 std::string(dbStateVarsLabels[i]));

  // Test _dbToStateVars()
  const int numLocs = _data->numLocs;
  scalar_array dbValues(numDBStateVars);

  const int stateVarsSize = _data->numVarsQuadPt;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBStateVars; ++i)
      dbValues[i] = _data->dbStateVars[iLoc*numDBStateVars+i];

    _material->_dbToStateVars(&stateVars[0], dbValues);
    
    const PylithScalar* const stateVarsE = 
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsE) ||
		    (0 == stateVarsSize && 0 == stateVarsE) );
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
pylith::materials::TestMaterial::testNonDimStateVars(void)
{ // testNonDimStateVars
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsQuadPt;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&stateVars[0], &_data->stateVars[iLoc*stateVarsSize],
	   stateVarsSize*sizeof(PylithScalar));
    _material->_nondimStateVars(&stateVars[0], stateVars.size());
    
    const PylithScalar* const stateVarsNondimE =
      (stateVarsSize > 0) ? &_data->stateVarsNondim[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsNondimE) ||
		    (0 == stateVarsSize && 0 == stateVarsNondimE) );

    const PylithScalar tolerance = 1.0e-06;
    for (int i=0; i < stateVarsSize; ++i) {
      if (fabs(stateVarsNondimE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsNondimE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsNondimE[i], stateVars[i],
				     tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testNonDimStateVars

// ----------------------------------------------------------------------
// Test _dimStateVars().
void
pylith::materials::TestMaterial::testDimStateVars(void)
{ // testDimStateVars
  PYLITH_METHOD_BEGIN;
 
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsQuadPt;
  scalar_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&stateVars[0], &_data->stateVarsNondim[iLoc*stateVarsSize],
	   stateVarsSize*sizeof(PylithScalar));
    _material->_dimStateVars(&stateVars[0], stateVars.size());
    
    const PylithScalar* const stateVarsE =
      (stateVarsSize > 0) ? &_data->stateVars[iLoc*stateVarsSize] : 0;
    CPPUNIT_ASSERT( (0 < stateVarsSize && stateVarsE) ||
		    (0 == stateVarsSize && 0 == stateVarsE) );

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


// End of file 
