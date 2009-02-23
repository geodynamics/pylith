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
//#include "data/ElasticIsotropic3DData.hh" // USES ElasticIsotropic3DData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

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

#if 0
// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestElasticMaterial::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  Quadrature<topology::Mesh> quadrature;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &quadrature, &material, &data);

  CPPUNIT_ASSERT(false);
} // testInitialize

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticMaterial::testCalcDensity(void)
{ // testCalcDensity
  topology::Mesh mesh;
  Quadrature<topology::Mesh> quadrature;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &quadrature, &material, &data);

  // Get cells associated with material
  const ALE::Obj<SieveMesh::label_sequence>& cells = mesh->heightStratum(0);
  SieveMesh::label_sequence::iterator c_iter = cells->begin();
  material.retrievePropsAndVars(*c_iter);
  const double_array& density = material.calcDensity();

  CPPUNIT_ASSERT_EQUAL(data._numQuadPts, density.size());

  const double tolerance = 1.0e-06;
  const double* densityE = data.density;
  CPPUNIT_ASSERT(0 != densityE);
  const double size = data._numQuadPts;
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
#endif    

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

  const int numLocs = _data->numLocs;
  const int numPropsQuadPt = _data->numPropsQuadPt;
  const int numVarsQuadPt = _data->numVarsQuadPt;
  
  double density = 0;
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &_data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &_data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));

    _matElastic->_calcDensity(&density, 
			      &properties[0], properties.size(),
			      &stateVars[0], stateVars.size());
    
    const double densityE = data->density[iLoc];
    
    const double tolerance = 1.0e-06;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density/densityE, tolerance);
  } // for
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcStress()
void
pylith::materials::TestElasticMaterial::test_calcStress(void)
{ // _testCalcStress
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const bool computeStateVars = true;

  const int numLocs = _data->numLocs;
  const int numPropsQuadPt = _data->numPropsQuadPt;
  const int numVarsQuadPt = _data->numVarsQuadPt;
  const int tensorSize = _matElastic->_tensorSize;
  
  double_array stress(tensorSize);
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);
  double_array strain(tensorSize);
  double_array initialStress(tensorSize);
  double_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(double));

    _matElastic->_calcStress(&stress[0], stress.size(),
			     &properties[0], properties.size(),
			     &stateVars[0], stateVars.size(),
			     &strain[0], strain.size(),
			     &initialStress[0], initialStress.size(),
			     &initialStrain[0], initialStrain.size(),
			     computeStateVars);
  
    const double* stressE = &data->stress[iLoc*tensorSize];
    CPPUNIT_ASSERT(0 != stressE);

    const double tolerance = 1.0e-06;
    for (int i=0; i < tensorSize; ++i)
      if (fabs(stressE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[i],
				     tolerance);
  } // for
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
  int tensorSize = 0;
  switch(data->dimension)
    { // switch
    case 1 :
      numConsts = 1;
      tensorSize = 1;
      break;
    case 2 :
      numConsts = 6;
      tensorSize = 3;
      break;
    case 3 :
      numConsts = 21;
      tensorSize = 6;
      break;
    } // switch
  const int numLocs = _data->numLocs;
  const int numPropsQuadPt = _data->numPropsQuadPt;
  const int numVarsQuadPt = _data->numVarsQuadPt;
  
  double_array elasticConsts(numConsts);
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);
  double_array strain(tensorSize);
  double_array initialStress(tensorSize);
  double_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(double));

    _matElastic->_calcElasticConsts(&elasticConsts[0], elasticConsts.size(),
				    &properties[0], properties.size(),
				    &stateVars[0], stateVars.size(),
				    &strain[0], strain.size(),
				    &initialStress[0], initialStress.size(),
				    &initialStrain[0], initialStrain.size());

    const double* elasticConstsE = &data->elasticConsts[iLoc*numConsts];
    CPPUNIT_ASSERT(0 != elasticConstsE);
    
    const double tolerance = 1.0e-06;
    for (int i=0; i < numConsts; ++i)
      if (fabs(elasticConstsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i],
				     tolerance);
  } // for
} // _testCalcElasticConsts

// ----------------------------------------------------------------------
// Test _updateStateVars()
void
pylith::materials::TestElasticMaterial::test_updateStateVars(void)
{ // test_updateStateVars
  CPPUNIT_ASSERT(0 != _matElastic);
  CPPUNIT_ASSERT(0 != _dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const bool computeStateVars = true;

  const int numLocs = _data->numLocs;
  const int numPropsQuadPt = _data->numPropsQuadPt;
  const int numVarsQuadPt = _data->numVarsQuadPt;
  const int tensorSize = _matElastic->_tensorSize;
  
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);
  double_array strain(tensorSize);
  double_array initialStress(tensorSize);
  double_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(double));

    _matElastic->_updateStateVars(&stateVars[0], stateVars.size(),
				  &properties[0], properties.size(),
				  &strain[0], strain.size(),
				  &initialStress[0], initialStress.size(),
				  &initialStrain[0], initialStrain.size());
    
    const double* stateVarsE = 
      (numVarsQuadPt > 0) ? &data->stateVarsUpdated[iLoc*numVarsQuadPt] : 0;
    CPPUNIT_ASSERT( (0 < numVarsQuadPt && 0 != stateVarsE) ||
		    (0 == numVarsQuadPt && 0 == stateVarsE) );

    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsQuadPt; ++i)
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i],
				     tolerance);
  } // for
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
    _matElastic->_stableTimeStepImplicit(data->properties, data->numPropsQuadPt,
					 data->stateVars, data->numVarsQuadPt);

  const double dtE = data->dtStableImplicit;

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);
} // _testStableTimeStepImplicit

// ----------------------------------------------------------------------
// Setup nondimensionalization.
void
pylith::materials::TestElasticMaterial::setupNormalizer(void)
{ // setupNormalizer
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.timeScale(_data->timeScale);
  normalizer.densityScale(_data->densityScale);
  _material->normalizer(normalizer);
  _matElastic->normalizer(normalizer);
} // setupNormalizer


// End of file 
