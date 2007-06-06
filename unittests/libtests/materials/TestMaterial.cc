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
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <petscmesh.h> // USES PETSc Mesh

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
    ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
					 const_cast<int*>(cells), numVertices,
					   interpolate, numCorners);
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
  const int cellDim = 1;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 1;
  const double verticesRef[] = { -1.0, 1.0, 0.0 };
  const double basis[] = { 0.455, 0.667, -0.122, -0.122, 0.667, 0.455 };
  const double basisDeriv[] = { -1.077, 1.155, -0.077, 0.077, -1.155, 1.077 };
  const double quadPtsRef[] = { -0.577350269, 0.577350269 };
  const double quadWts[] = { 1.0, 1.0  };
  quadrature.initialize(verticesRef, basis, basisDeriv, quadPtsRef, quadWts,
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

  Mesh::label_sequence::iterator cellIter=cells->begin();
  const double tolerance = 1.0e-06;

  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  const real_section_type::value_type* densityCell = 
    parameterDensity->restrictPoint(*cellIter);
  CPPUNIT_ASSERT(0 != densityCell);
  for (int i=0; i < numQuadPts; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, densityCell[i]/densityE[i], tolerance);
  
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  const real_section_type::value_type* muCell = 
    parameterMu->restrictPoint(*cellIter);
  CPPUNIT_ASSERT(0 != muCell);
  for (int i=0; i < numQuadPts; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, muCell[i]/muE[i], tolerance);
  
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  const real_section_type::value_type* lambdaCell = 
    parameterLambda->restrictPoint(*cellIter);
  CPPUNIT_ASSERT(0 != lambdaCell);
  for (int i=0; i < numQuadPts; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, lambdaCell[i]/lambdaE[i], tolerance);
} // testInitialize

// ----------------------------------------------------------------------
// Test DBToParameters
void
pylith::materials::TestMaterial::_testDBToParameters(Material* material,
						     const MaterialData& data) const
{ // _testDBToParameters
  CPPUNIT_ASSERT(0 != material);
  
  const int numLocs = data.numLocs;
  const int numDBValues = data.numDBValues;
  double_array dbData(numDBValues);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    for (int i=0; i < numDBValues; ++i)
      dbData[i] = data.dbData[iLoc*numDBValues+i];

    const int numParameters = data.numParameters;
    int numParamEntries = 0;

    std::vector<double_array> parameterData(numParameters);
    for (int iParam=0; iParam < numParameters; ++iParam) {
      parameterData[iParam].resize(data.numParamValues[iParam]);
      numParamEntries += data.numParamValues[iParam];
    } // for

    double* const parameterDataE = &data.parameterData[iLoc*numParamEntries];
    material->_dbToParameters(&parameterData, dbData);

    const double tolerance = 1.0e-06;
    for (int iParam=0, i=0; iParam < numParameters; ++iParam) {
      const int numParamValues = data.numParamValues[iParam];
      CPPUNIT_ASSERT_EQUAL(numParamValues, int(parameterData[iParam].size()));
      for (int iValue=0; iValue < numParamValues; ++iValue)
	if (fabs(parameterDataE[i]) > tolerance)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				       parameterData[iParam][iValue]/parameterDataE[i++],
				       tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(parameterDataE[i++], parameterData[iParam][iValue],
				       tolerance);
    } // for
  } // for
} // _testDBToParameters

// ----------------------------------------------------------------------
// Test _dbValues() and _numDBValues()
void
pylith::materials::TestMaterial::_testDBValues(Material* material,
					       const MaterialData& data) const
{ // _testDBValues
  CPPUNIT_ASSERT(0 != material);

  const int numDBValues = data.numDBValues;

  CPPUNIT_ASSERT(numDBValues == material->_numDBValues());
  char** const dbValuesE = data.dbValues;
  const char** dbValues = material->_dbValues();
  for (int i=0; i < numDBValues; ++i)
    CPPUNIT_ASSERT(0 == strcmp(dbValuesE[i], dbValues[i]));
} // _testDBValues

// ----------------------------------------------------------------------
// Test _numParameters() and _parameterNames()
void
pylith::materials::TestMaterial::_testParameters(Material* material,
						 const MaterialData& data) const
{ // _testParameters
  CPPUNIT_ASSERT(0 != material);

  const int numParameters = data.numParameters;

  int_array numParamValues;
  material->_numParamValues(&numParamValues);

  CPPUNIT_ASSERT_EQUAL(numParameters, int(numParamValues.size()));
  char** const namesE = data.parameterNames;
  const char** names = material->_parameterNames();
  for (int i=0; i < numParameters; ++i) {
    CPPUNIT_ASSERT_EQUAL(data.numParamValues[i], numParamValues[i]);
    CPPUNIT_ASSERT(0 == strcmp(namesE[i], names[i]));
  } // for
} // _testParameters


// End of file 
