// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMaterial.hh" // Implementation of class methods

#include "data/MaterialData.hh" // USES MaterialData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "pylith/utils/sievetypes.hh" // USES Mesh

#include <string.h> // USES strcmp()
#include <math.h> // USES assert()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaterial );

// ----------------------------------------------------------------------
// Test db()
void
pylith::materials::TestMaterial::testDB(void)
{ // testDB
  const char* label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label);
  
  ElasticIsotropic3D material;
  material.db(&db);
  
  CPPUNIT_ASSERT(0 != material._db);
  CPPUNIT_ASSERT(0 == strcmp(label, material._db->label()));
} // testDB

// ----------------------------------------------------------------------
// Test initialStateDB()
void
pylith::materials::TestMaterial::testInitialStateDB(void)
{ // testInitialStateDB
  const char* label = "my_database";
  spatialdata::spatialdb::SimpleDB initialStateDB;
  initialStateDB.label(label);
  
  ElasticIsotropic3D material;
  material.initialStateDB(&initialStateDB);
  
  CPPUNIT_ASSERT(0 != material._initialStateDB);
  CPPUNIT_ASSERT(0 == strcmp(label, material._initialStateDB->label()));
} // testInitialStateDB
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
  const char* label = "the_database";
  ElasticIsotropic3D material;
  material.label(label);
  
  CPPUNIT_ASSERT(0 == strcmp(label, material.label().c_str()));
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
// Test initialize()
void
pylith::materials::TestMaterial::testInitialize(void)
{ // testInitialize
  ALE::Obj<Mesh> mesh;
  const int materialID = 24;
  { // create mesh
    const int cellDim = 1;
    const int numCorners = 3;
    const int spaceDim = 1;
    const int numVertices = 3;
    const int numCells = 1;
    const double vertCoords[] = { -1.0, 0.0, 1.0};
    const int cells[] = { 0, 1, 2 };
    CPPUNIT_ASSERT(0 != vertCoords);
    CPPUNIT_ASSERT(0 != cells);

    mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

    const bool interpolate = false;
    ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

    ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
                                                const_cast<int*>(cells), numVertices,
                                                interpolate, numCorners);
    std::map<Mesh::point_type,Mesh::point_type> renumbering;
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);

  } // create mesh

  { // set material ids
    const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
    const ALE::Obj<Mesh::label_type>& labelMaterials = mesh->createLabel("material-id");
    int i = 0;
    for(Mesh::label_sequence::iterator e_iter = cells->begin();
	e_iter != cells->end();
	++e_iter)
      mesh->setValue(labelMaterials, *e_iter, materialID);
  } // set material ids

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(1);
  cs.initialize();

  feassemble::Quadrature1D quadrature;
  feassemble::GeometryLine1D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 1;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 1;
  const double basis[] = { 0.455, 0.667, -0.122, -0.122, 0.667, 0.455 };
  const double basisDeriv[] = { -1.077, 1.155, -0.077, 0.077, -1.155, 1.077 };
  const double quadPtsRef[] = { -0.577350269, 0.577350269 };
  const double quadWts[] = { 1.0, 1.0  };
  quadrature.initialize(basis, basisDeriv, quadPtsRef, quadWts,
			cellDim, numCorners, numQuadPts, spaceDim);

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii iohandler;
  iohandler.filename("data/matinitialize.spatialdb");
  db.ioHandler(&iohandler);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  ElasticIsotropic3D material;
  material.db(&db);
  material.id(materialID);
  material.label("my_material");
  material.initialize(mesh, &cs, &quadrature);

  const double densityA = 2000.0;
  const double vsA = 100.0;
  const double vpA = 180.0;
  const double muA = vsA*vsA*densityA;
  const double lambdaA = vpA*vpA*densityA - 2.0*muA;
  const double densityB = 3000.0;
  const double vsB = 200.0;
  const double vpB = 400.0;
  const double muB = vsB*vsB*densityB;
  const double lambdaB = vpB*vpB*densityB - 2.0*muB;
  const double densityE[] = { densityA, densityB };
  const double muE[] = { muA, muB };
  const double lambdaE[] = { lambdaA, lambdaB };

  // Get cells associated with material
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  Mesh::label_sequence::iterator c_iter = cells->begin();
  const double tolerance = 1.0e-06;

  const real_section_type::value_type* paramsCell =
    material._properties->restrictPoint(*c_iter);
  CPPUNIT_ASSERT(0 != paramsCell);

  const int pidDensity = 0;
  const int pidMu = 1;
  const int pidLambda = 2;

  // density
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._totalPropsQuadPt + pidDensity;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, paramsCell[index]/densityE[i], tolerance);
  } // for
  
  // mu
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._totalPropsQuadPt + pidMu;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, paramsCell[index]/muE[i], tolerance);
  } // for
  
  // lambda
  for (int i=0; i < numQuadPts; ++i) {
    const int index = i*material._totalPropsQuadPt + pidLambda;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, paramsCell[index]/lambdaE[i], tolerance);
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
// Test DBToProperties
void
pylith::materials::TestMaterial::testDBToProperties(void)
{ // testDBToProperties
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);
  
  const int numLocs = _data->numLocs;
  const int numDBValues = _data->numDBValues;
  double_array dbData(numDBValues);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBValues; ++i)
      dbData[i] = _data->dbData[iLoc*numDBValues+i];

    const int numProperties = _data->numParameters;
    int numParamEntries = 0;
    for (int iParam=0; iParam < numProperties; ++iParam)
      numParamEntries += _data->numParamValues[iParam];

    double_array parameterData(numParamEntries);

    double* const parameterDataE = &_data->parameterData[iLoc*numParamEntries];
    _material->_dbToProperties(&parameterData[0], dbData);

    const double tolerance = 1.0e-06;
    for (int i=0; i < numParamEntries; ++i) {
      if (fabs(parameterDataE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     parameterData[i]/parameterDataE[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(parameterDataE[i], parameterData[i],
				     tolerance);
    } // for
  } // for
} // testDBToProperties

// ----------------------------------------------------------------------
// Test _dbValues and _numDBValues.
void
pylith::materials::TestMaterial::testDBValues(void)
{ // testDBValues
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);

  const int numDBValues = _data->numDBValues;
  CPPUNIT_ASSERT_EQUAL(numDBValues, _material->_numDBValues);

  char** const dbValuesE = _data->dbValues;
  for (int i=0; i < numDBValues; ++i)
    CPPUNIT_ASSERT(0 == strcmp(dbValuesE[i], _material->_dbValues[i]));
} // testDBValues

// ----------------------------------------------------------------------
// Test _numProperties.
void
pylith::materials::TestMaterial::testProperties(void)
{ // testProperties
  CPPUNIT_ASSERT(0 != _material);
  CPPUNIT_ASSERT(0 != _data);

  const int numProperties = _data->numParameters;
  CPPUNIT_ASSERT_EQUAL(numProperties, _material->_numProperties);
} // testProperties


// End of file 
