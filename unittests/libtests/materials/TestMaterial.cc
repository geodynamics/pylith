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

#include "TestMaterial.hh" // Implementation of class methods

#include "data/MaterialData.hh" // USES MaterialData

#include "pylith/utils/array.hh" // USES double_array

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcmp()
#include <cassert> // USES assert()
#include <cmath> // USES sqrt()

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaterial );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Test id()
void
pylith::materials::TestMaterial::testID(void)
{ // testID
  const int id = 346;
  ElasticIsotropic3D material;
  material.id(id);
  
  CPPUNIT_ASSERT(id == material.id());
} // testID

// ----------------------------------------------------------------------
// Test label()
void
pylith::materials::TestMaterial::testLabel(void)
{ // testLabel
  const std::string& label = "the database";
  ElasticIsotropic3D material;
  material.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(material.label()));
} // testLabel
    
// ----------------------------------------------------------------------
// Test timestep()
void
pylith::materials::TestMaterial::testTimeStep(void) 
{ // testTimeStep
  const double dt = 2.0;
  ElasticIsotropic3D material;
  material.timeStep(dt);
  
  CPPUNIT_ASSERT_EQUAL(dt, material.timeStep());
} // testTimeStep

// ----------------------------------------------------------------------
// Test dbProperties()
void
pylith::materials::TestMaterial::testDBProperties(void)
{ // testDBProperties
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticIsotropic3D material;
  material.dbProperties(&db);
  
  CPPUNIT_ASSERT(0 != material._dbProperties);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbProperties->label()));
} // testDBProperties

// ----------------------------------------------------------------------
// Test dbStateVars()
void
pylith::materials::TestMaterial::testDBStateVars(void)
{ // testDBStateVars
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticIsotropic3D material;
  material.dbInitialState(&db);
  
  CPPUNIT_ASSERT(0 != material._dbInitialState);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialState->label()));
} // testDBStateVars

// ----------------------------------------------------------------------
// Test normalizer()
void
pylith::materials::TestMaterial::testNormalizer(void)
{ // testNormalizer
  spatialdata::units::Nondimensional normalizer;
  const double lengthScale = 2.0;
  normalizer.lengthScale(lengthScale);

  ElasticIsotropic3D material;
  material.normalizer(normalizer);
  
  CPPUNIT_ASSERT(0 != material._normalizer);
  CPPUNIT_ASSERT_EQUAL(lengthScale, material._normalizer->lengthScale());
} // testNormalizer

// ----------------------------------------------------------------------
// Test needNewJacobian()
void
pylith::materials::TestMaterial::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  ElasticIsotropic3D material;

  bool flag = false;
  material._needNewJacobian = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.needNewJacobian());

  flag = true;
  material._needNewJacobian = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.needNewJacobian());
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test isJacobianSymmetric()
void
pylith::materials::TestMaterial::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  ElasticIsotropic3D material;

  CPPUNIT_ASSERT_EQUAL(true, material.isJacobianSymmetric());

  bool flag = false;
  material._isJacobianSymmetric = flag;
  CPPUNIT_ASSERT_EQUAL(flag, material.isJacobianSymmetric());
} // testIsJacobianSymmetric

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestMaterial::testInitialize(void)
{ // testInitialize
  // Setup mesh
  topology::Mesh mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/line3.mesh");
  iohandler.read(&mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);

  // Setup quadrature
  feassemble::Quadrature<topology::Mesh> quadrature;
  feassemble::GeometryLine1D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 1;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 1;
  const double basis[] = { 0.455, -0.122, 0.667, -0.122, 0.455, 0.667 };
  const double basisDeriv[] = { 
    -1.07735027e+00,
    -7.73502692e-02,
    1.15470054e+00,
    7.73502692e-02,
    1.07735027e+00,
    -1.15470054e+00,
  };
  const double quadPtsRef[] = { -0.577350269, 0.577350269 };
  const double quadWts[] = { 1.0, 1.0  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);


  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);

  // Compute geometry for cells
  quadrature.initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  quadrature.computeGeometry(mesh, cells);
#endif

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename("data/matinitialize.spatialdb");
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  spatialdata::units::Nondimensional normalizer;

  ElasticStrain1D material;
  material.dbProperties(&db);
  material.id(materialId);
  material.label("my_material");
  material.normalizer(normalizer);
  material.initialize(mesh, &quadrature);

  const double densityA = 2500.0;
  const double vsA = 3000.0;
  const double vpA = vsA*sqrt(3.0);
  const double muA = vsA*vsA*densityA;
  const double lambdaA = vpA*vpA*densityA - 2.0*muA;
  const double densityB = 2000.0;
  const double vsB = 1200.0;
  const double vpB = vsB*sqrt(3.0);
  const double muB = vsB*vsB*densityB;
  const double lambdaB = vpB*vpB*densityB - 2.0*muB;
  const double densityE[] = { densityA, densityB };
  const double muE[] = { muA, muB };
  const double lambdaE[] = { lambdaA, lambdaB };

  SieveMesh::point_type cell = *cells->begin();
  const double tolerance = 1.0e-06;

  CPPUNIT_ASSERT(0 != material._properties);
  const Obj<RealSection>& propertiesSection = material._properties->section();
  CPPUNIT_ASSERT(!propertiesSection.isNull());
  const double* propertiesCell = propertiesSection->restrictPoint(cell);
  CPPUNIT_ASSERT(0 != propertiesCell);

  const int p_density = 0;
  const int p_mu = 1;
  const int p_lambda = 2;

  // density
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_density;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesCell[index]/densityE[i],
				 tolerance);
  } // for
  
  // mu
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_mu;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesCell[index]/muE[i], tolerance);
  } // for
  
  // lambda
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._numPropsQuadPt + p_lambda;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, propertiesCell[index]/lambdaE[i], 
				 tolerance);
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaterial::setUp(void)
{ // setUp
  _material = 0;
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::materials::TestMaterial::tearDown(void)
{ // tearDown
  delete _material; _material = 0;
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test dimension()
void
pylith::materials::TestMaterial::testDimension(void)
{ // testDimension
  CPPUNIT_ASSERT_EQUAL(_data->dimension, _material->dimension());
} // testDimension

// ----------------------------------------------------------------------
// Test dimension()
void
pylith::materials::TestMaterial::testTensorSize(void)
{ // testTensorSize
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
} // testTensorSize

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::materials::TestMaterial::testDBToProperties(void)
{ // testDBToProperties
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
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
  double_array dbValues(numDBProperties);

  const int propertiesSize = _data->numPropsQuadPt;
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBProperties; ++i)
      dbValues[i] = _data->dbProperties[iLoc*numDBProperties+i];

    _material->_dbToProperties(&properties[0], dbValues);
    
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
pylith::materials::TestMaterial::testNonDimProperties(void)
{ // testNonDimProperties
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsQuadPt;
  double_array propertiesNondim(propertiesSize);
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&properties[0], &_data->properties[iLoc*propertiesSize],
	   propertiesSize*sizeof(double));
    _material->_nondimProperties(&properties[0], properties.size());
    
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
pylith::materials::TestMaterial::testDimProperties(void)
{ // testDimProperties
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int propertiesSize = _data->numPropsQuadPt;
  double_array properties(propertiesSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&properties[0], &_data->propertiesNondim[iLoc*propertiesSize], 
	   propertiesSize*sizeof(double));
    _material->_dimProperties(&properties[0], properties.size());
    
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
pylith::materials::TestMaterial::testDBToStateVars(void)
{ // testDBToStateVars
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
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
  double_array dbValues(numDBStateVars);

  const int stateVarsSize = _data->numVarsQuadPt;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBStateVars; ++i)
      dbValues[i] = _data->dbStateVars[iLoc*numDBStateVars+i];

    _material->_dbToStateVars(&stateVars[0], dbValues);
    
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
pylith::materials::TestMaterial::testNonDimStateVars(void)
{ // testNonDimStateVars
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsQuadPt;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&stateVars[0], &_data->stateVars[iLoc*stateVarsSize],
	   stateVarsSize*sizeof(double));
    _material->_nondimStateVars(&stateVars[0], stateVars.size());
    
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
pylith::materials::TestMaterial::testDimStateVars(void)
{ // testDimStateVars
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int stateVarsSize = _data->numVarsQuadPt;
  double_array stateVars(stateVarsSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {

    memcpy(&stateVars[0], &_data->stateVarsNondim[iLoc*stateVarsSize],
	   stateVarsSize*sizeof(double));
    _material->_dimStateVars(&stateVars[0], stateVars.size());
    
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


// End of file 
