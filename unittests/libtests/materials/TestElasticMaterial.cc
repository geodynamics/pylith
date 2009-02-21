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

#include "TestElasticMaterial.hh" // Implementation of class methods

#include "data/ElasticMaterialData.hh" // USES ElasticMaterialData
#include "data/ElasticIsotropic3DData.hh" // USES ElasticIsotropic3DData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticMaterial );

// ----------------------------------------------------------------------
// Test dbInitialStress()
void
pylith::materials::TestElasticMaterial::testDBInitialStress(void)
{ // testDBInitialStress
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticIsotropic3D material;
  material.dbInitialStress(&db);
  
  CPPUNIT_ASSERT(0 != material._dbInitialStress);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialStress->label()));
} // testDBInitialStress

// ----------------------------------------------------------------------
// Test dbInitialStrain()
void
pylith::materials::TestElasticMaterial::testDBInitialStrain(void)
{ // testDBInitialStrain
  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticIsotropic3D material;
  material.dbInitialStrain(&db);
  
  CPPUNIT_ASSERT(0 != material._dbInitialStrain);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialStrain->label()));
} // testDBInitialStrain

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestElasticMaterial::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(false);
} // testInitialize

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticMaterial::testCalcDensity(void)
{ // testCalcDensity
  ALE::Obj<Mesh> mesh;
  { // create mesh
    const int cellDim = 1;
    const int numCorners = 2;
    const int spaceDim = 1;
    const int numVertices = 2;
    const int numCells = 1;
    const double vertCoords[] = { -1.0, 1.0};
    const int cells[] = { 0, 1};
    CPPUNIT_ASSERT(0 != vertCoords);
    CPPUNIT_ASSERT(0 != cells);

    mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

    const bool interpolate = false;
    ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

    ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    std::map<Mesh::point_type,Mesh::point_type> renumbering;
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  const int numQuadPts = 2;
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;
  
  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._properties = new real_section_type(mesh->comm(), mesh->debug());
  material._properties->setChart(mesh->getSieve()->getChart());
  material._properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._properties);

  material._properties->updatePoint(*c_iter, data.parameterData);

  material.getPropertiesCell(*c_iter, numQuadPts);
  const double_array& density = material.calcDensity();

  const double tolerance = 1.0e-06;
  const double* densityE = data.density;
  CPPUNIT_ASSERT(0 != densityE);
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
} // testCalcDensity
    
// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticMaterial::testCalcStress(void)
{ // testCalcProperties
  ALE::Obj<Mesh> mesh;
  { // create mesh
    const int cellDim = 1;
    const int numCorners = 2;
    const int spaceDim = 1;
    const int numVertices = 2;
    const int numCells = 1;
    const double vertCoords[] = { -1.0, 1.0};
    const int cells[] = { 0, 1};
    CPPUNIT_ASSERT(0 != vertCoords);
    CPPUNIT_ASSERT(0 != cells);

    mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

    const bool interpolate = false;
    ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

    ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    std::map<Mesh::point_type,Mesh::point_type> renumbering;
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  const int numQuadPts = 2;
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._properties = new real_section_type(mesh->comm(), mesh->debug());
  material._properties->setChart(mesh->getSieve()->getChart());
  material._properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._properties);

  material._properties->updatePoint(*c_iter, data.parameterData);
  const int strainSize = material._tensorSize;
  double_array strain(data.strain, numQuadPts*strainSize);

  const int initialStateSize = material._tensorSize;
  material._initialStateSize = initialStateSize;
  const int initialStateFiberDim = numQuadPts * initialStateSize;
  material._initialState = new real_section_type(mesh->comm(), mesh->debug());
  material._initialState->setChart(mesh->getSieve()->getChart());
  material._initialState->setFiberDimension(cells, initialStateFiberDim);
  mesh->allocate(material._initialState);
  material._initialState->updatePoint(*c_iter, data.initialState);

  material.getPropertiesCell(*c_iter, numQuadPts);
  const double_array& stress = material.calcStress(strain);

  const double tolerance = 1.0e-06;
  const double* stressE = data.stress;
  CPPUNIT_ASSERT(0 != stressE);
  
  const int size = numQuadPts * strainSize;;
  CPPUNIT_ASSERT_EQUAL(size, int(stress.size()));
  for (int i=0; i < size; ++i)
    if (fabs(stressE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[i], tolerance);
} // testCalcStress
    
// ----------------------------------------------------------------------
// Test calcDerivElastic()
void
pylith::materials::TestElasticMaterial::testCalcDerivElastic(void)
{ // testCalcDerivElastic
  ALE::Obj<Mesh> mesh;
  { // create mesh
    const int cellDim = 1;
    const int numCorners = 2;
    const int spaceDim = 1;
    const int numVertices = 2;
    const int numCells = 1;
    const double vertCoords[] = { -1.0, 1.0};
    const int cells[] = { 0, 1};
    CPPUNIT_ASSERT(0 != vertCoords);
    CPPUNIT_ASSERT(0 != cells);

    mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

    const bool interpolate = false;
    ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

    ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    std::map<Mesh::point_type,Mesh::point_type> renumbering;
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  const int numQuadPts = 2;
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._properties = new real_section_type(mesh->comm(), mesh->debug());
  material._properties->setChart(mesh->getSieve()->getChart());
  material._properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._properties);

  material._properties->updatePoint(*c_iter, data.parameterData);
  const int strainSize = material._tensorSize;
  double_array strain(data.strain, numQuadPts*strainSize);

  const int initialStateSize = material._tensorSize;
  material._initialStateSize = initialStateSize;
  const int initialStateFiberDim = numQuadPts * initialStateSize;
  material._initialState = new real_section_type(mesh->comm(), mesh->debug());
  material._initialState->setChart(mesh->getSieve()->getChart());
  material._initialState->setFiberDimension(cells, initialStateFiberDim);
  mesh->allocate(material._initialState);

  material.getPropertiesCell(*c_iter, numQuadPts);
  const double_array& elasticConsts = material.calcDerivElastic(strain);

  const double tolerance = 1.0e-06;
  const double* elasticConstsE = data.elasticConsts;
  CPPUNIT_ASSERT(0 != elasticConstsE);

  const int size = elasticConsts.size();
  for (int i=0; i < size; ++i)
    if (fabs(elasticConstsE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				   elasticConsts[i]/elasticConstsE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], 
				   elasticConsts[i],
				   tolerance);
} // testCalcDerivElastic
    
// ----------------------------------------------------------------------
// Test updateStateVars()
void
pylith::materials::TestElasticMaterial::testUpdateStateVars(void)
{ // testUpdateStateVars
  CPPUNIT_ASSERT(false);
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcStableTimeStepImplicit()
void
pylith::materials::TestElasticMaterial::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  ALE::Obj<Mesh> mesh;
  { // create mesh
    const int cellDim = 1;
    const int numCorners = 2;
    const int spaceDim = 1;
    const int numVertices = 2;
    const int numCells = 1;
    const double vertCoords[] = { -1.0, 1.0};
    const int cells[] = { 0, 1};
    CPPUNIT_ASSERT(0 != vertCoords);
    CPPUNIT_ASSERT(0 != cells);

    mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

    const bool interpolate = false;
    ALE::Obj<ALE::Mesh::sieve_type> s = 
      new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

    ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    std::map<Mesh::point_type,Mesh::point_type> renumbering;
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  const int numQuadPts = 2;
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  const int initialStateSize = material._tensorSize;
  material._initialStateSize = initialStateSize;
  const int initialStateFiberDim = numQuadPts * initialStateSize;
  material._initialState = new real_section_type(mesh->comm(), mesh->debug());
  material._initialState->setChart(mesh->getSieve()->getChart());
  material._initialState->setFiberDimension(cells, initialStateFiberDim);
  mesh->allocate(material._initialState);
  
  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._properties = new real_section_type(mesh->comm(), mesh->debug());
  material._properties->setChart(mesh->getSieve()->getChart());
  material._properties->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._properties);

  material._properties->updatePoint(*c_iter, data.parameterData);

  material.getPropertiesCell(*c_iter, numQuadPts);
  const double dt = material.stableTimeStepImplicit();

  const double tolerance = 1.0e-06;
  const double dtE = data.dtStableImplicit;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);
} // testStableTimeStepImplicit
    
// ----------------------------------------------------------------------
// Test useElasticBehavior()
void
pylith::materials::TestElasticMaterial::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  ElasticIsotropic3D material;

  bool flag = false;
  material.useElasticBehavior(flag);
  CPPUNIT_ASSERT_EQUAL(flag, material._useElasticBehavior);

  bool flag = true;
  material.useElasticBehavior(flag);
  CPPUNIT_ASSERT_EQUAL(flag, material._useElasticBehavior);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticMaterial::setUp(void)
{ // setUp
  TestMaterial::setUp();
  _matElastic = 0;
  _dataElastic = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::materials::TestElasticMaterial::tearDown(void)
{ // tearDown
  TestMaterial::tearDown();
  delete _matElastic; _matElastic = 0;
  delete _dataElastic; _dataElastic = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test _calcDensity()
void
pylith::materials::TestElasticMaterial::test_calcDensity(void)
{ // _testCalcDensity
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  double_array density(1);
  _matElastic->_calcDensity(&density[0], data->parameterData, data->numParamsQuadPt);

  const double* densityE = data->density;
  CPPUNIT_ASSERT(0 != densityE);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[0]/densityE[0], tolerance);
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcStress()
void
pylith::materials::TestElasticMaterial::test_calcStress(void)
{ // _testCalcElasticConsts
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const bool computeStateVars = true;

  const int stressSize = _matElastic->_tensorSize;
  double_array stress(stressSize);
  _matElastic->_calcStress(&stress[0], stress.size(),
			 data->parameterData, data->numParamsQuadPt,
			   data->strain, stressSize, 
			   data->initialState, data->numInitialStateValues,
			   computeStateVars);
  
  const double* stressE = data->stress;
  CPPUNIT_ASSERT(0 != stressE);

  const double tolerance = 1.0e-06;
  for (int i=0; i < stressSize; ++i)
    if (fabs(stressE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[i],
				   tolerance);
} // _testCalcStress

// ----------------------------------------------------------------------
// Test _calcElasticConsts()
void
pylith::materials::TestElasticMaterial::test_calcElasticConsts(void)
{ // _testCalcElasticConsts
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  int numConsts = 0;
  int strainSize = 0;
  switch(data->dimension)
    { // switch
    case 1 :
      numConsts = 1;
      strainSize = 1;
      break;
    case 2 :
      numConsts = 6;
      strainSize = 3;
      break;
    case 3 :
      numConsts = 21;
      strainSize = 6;
      break;
    } // switch
  double_array elasticConsts(numConsts);
  _matElastic->_calcElasticConsts(&elasticConsts[0], numConsts,
				data->parameterData, data->numParamsQuadPt,
				data->strain, strainSize,
				data->initialState,
				data->numInitialStateValues);

  const double* elasticConstsE = data->elasticConsts;
  CPPUNIT_ASSERT(0 != elasticConstsE);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numConsts; ++i)
    if (fabs(elasticConstsE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i],
				   tolerance);
} // _testCalcElasticConsts

// ----------------------------------------------------------------------
// Test _updateStateVars()
void
pylith::materials::TestElasticMaterial::test_updateStateVars(void)
{ // test_updateStateVars
  CPPUNIT_ASSERT(false);
} // test_updateStateVars

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestElasticMaterial::test_stableTimeStepImplicit(void)
{ // _testCalcDensity
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const double dt =
  _matElastic->_stableTimeStepImplicit(data->parameterData, 
				       data->numParamsQuadPt);

  const double dtE = data->dtStableImplicit;

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);
} // _testStableTimeStepImplicit


// End of file 
