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

#include "TestElasticMaterial.hh" // Implementation of class methods

#include "data/ElasticMaterialData.hh" // USES ElasticMaterialData
#include "data/ElasticStrain1DData.hh" // USES ElasticStrain1DData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticMaterial );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  const PylithScalar tolerance = 1.0e-06;
  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Test initialStress field
  CPPUNIT_ASSERT(0 != material._initialFields);
  const ALE::Obj<RealSection>& stressSection = 
    material._initialFields->get("initial stress").section();
  CPPUNIT_ASSERT(!stressSection.isNull());
  int fiberDim = numQuadPts * tensorSize;
  CPPUNIT_ASSERT_EQUAL(fiberDim, stressSection->getFiberDimension(cell));
  const PylithScalar* initialStress = stressSection->restrictPoint(cell);
  CPPUNIT_ASSERT(0 != initialStress);
  const PylithScalar* initialStressE = data.initialStress;
  CPPUNIT_ASSERT(0 != initialStressE);
  for (int i=0; i < fiberDim; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStress[i]/initialStressE[i],
				 tolerance);

  // Test initialStrain field
  CPPUNIT_ASSERT(0 != material._initialFields);
  const ALE::Obj<RealSection>& strainSection = 
    material._initialFields->get("initial strain").section();
  CPPUNIT_ASSERT(!strainSection.isNull());
  fiberDim = numQuadPts * tensorSize;
  CPPUNIT_ASSERT_EQUAL(fiberDim, strainSection->getFiberDimension(cell));
  const PylithScalar* initialStrain = strainSection->restrictPoint(cell);
  CPPUNIT_ASSERT(0 != initialStrain);
  const PylithScalar* initialStrainE = data.initialStrain;
  CPPUNIT_ASSERT(0 != initialStrainE);
  for (int i=0; i < fiberDim; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStrain[i]/initialStrainE[i],
				 tolerance);

  // Test cell arrays
  size_t size = data.numLocs*data.numPropsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, material._propertiesCell.size());

  size = data.numLocs*data.numVarsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, material._stateVarsCell.size());

  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, material._initialStressCell.size());

  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, material._initialStrainCell.size());

  size = data.numLocs;
  CPPUNIT_ASSERT_EQUAL(size, material._densityCell.size());

  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, material._stressCell.size());

  int numElasticConsts = 0;
  switch (data.dimension)
    { // switch
    case 1 :
      numElasticConsts = 1;
      break;
    case 2 :
      numElasticConsts = 6;
      break;
    case 3 :
      numElasticConsts = 21;
      break;
    default :
      assert(0);
    } // switch
  size = data.numLocs*numElasticConsts;
  CPPUNIT_ASSERT_EQUAL(size, material._elasticConstsCell.size());
} // testInitialize

// ----------------------------------------------------------------------
// Test retrievePropsAndVars().
void
pylith::materials::TestElasticMaterial::testRetrievePropsAndVars(void)
{ // testRetrievePropsAndVars
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  material.retrievePropsAndVars(cell);

  const PylithScalar tolerance = 1.0e-06;
  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;
  const int numVarsQuadPt = data.numVarsQuadPt;

  // Test cell arrays
  const PylithScalar* propertiesE = data.properties;
  CPPUNIT_ASSERT(0 != propertiesE);
  const scalar_array& properties = material._propertiesCell;
  size_t size = data.numLocs*data.numPropsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, properties.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, properties[i]/propertiesE[i],
				 tolerance);

  const PylithScalar* stateVarsE = data.stateVars;
  CPPUNIT_ASSERT( (0 < numVarsQuadPt && 0 != stateVarsE) ||
		  (0 == numVarsQuadPt && 0 == stateVarsE) );
  const scalar_array& stateVars = material._stateVarsCell;
  size = data.numLocs*numVarsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, stateVars.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i],
				 tolerance);

  const PylithScalar* initialStressE = data.initialStress;
  CPPUNIT_ASSERT(0 != initialStressE);
  const scalar_array& initialStress = material._initialStressCell;
  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, initialStress.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStress[i]/initialStressE[i],
				 tolerance);

  const PylithScalar* initialStrainE = data.initialStrain;
  CPPUNIT_ASSERT(0 != initialStrainE);
  const scalar_array& initialStrain = material._initialStrainCell;
  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, initialStrain.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStrain[i]/initialStrainE[i],
				 tolerance);
} // testRetrievePropsAndVars

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticMaterial::testCalcDensity(void)
{ // testCalcDensity
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  material.retrievePropsAndVars(cell);
  const scalar_array& density = material.calcDensity();

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  const PylithScalar* densityE = data.density;
  CPPUNIT_ASSERT(0 != densityE);
  const size_t size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, density.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
} // testCalcDensity
    
// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticMaterial::testCalcStress(void)
{ // testCalcStress
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Setup total strain
  scalar_array strain(data.strain, numQuadPts*tensorSize);

  material.retrievePropsAndVars(cell);
  const scalar_array& stress = material.calcStress(strain);

  const PylithScalar* stressE = data.stress;
  CPPUNIT_ASSERT(0 != stressE);
  const size_t size = numQuadPts * tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, stress.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], tolerance);
} // testCalcStress
    
// ----------------------------------------------------------------------
// Test calcDerivElastic()
void
pylith::materials::TestElasticMaterial::testCalcDerivElastic(void)
{ // testCalcDerivElastic
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Setup total strain
  scalar_array strain(data.strain, numQuadPts*tensorSize);

  material.retrievePropsAndVars(cell);
  const scalar_array& elasticConsts = material.calcDerivElastic(strain);

  int numElasticConsts = 0;
  switch (data.dimension)
    { // switch
    case 1 :
      numElasticConsts = 1;
      break;
    case 2 :
      numElasticConsts = 9;
      break;
    case 3 :
      numElasticConsts = 36;
      break;
    default :
      assert(0);
    } // switch

  const PylithScalar* elasticConstsE = data.elasticConsts;
  CPPUNIT_ASSERT(0 != elasticConstsE);
  const size_t size = numQuadPts * numElasticConsts;
  CPPUNIT_ASSERT_EQUAL(size, elasticConsts.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i],
				 tolerance);
} // testCalcDerivElastic
    
// ----------------------------------------------------------------------
// Test updateStateVars()
void
pylith::materials::TestElasticMaterial::testUpdateStateVars(void)
{ // testUpdateStateVars
  std::cout << "\n\nWARNING!! WARNING!! WARNING!!\n"
    "Need to implement using material with state variables.\n\n";
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcStableTimeStepImplicit()
void
pylith::materials::TestElasticMaterial::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  topology::Mesh mesh;
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  SieveMesh::point_type cell = *cells->begin();

  material.retrievePropsAndVars(cell);
  const PylithScalar dt = material.stableTimeStepImplicit(mesh);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dtE = data.dtStableImplicit;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);
} // testStableTimeStepImplicit

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

  const int numLocs = data->numLocs;
  const int numPropsQuadPt = data->numPropsQuadPt;
  const int numVarsQuadPt = data->numVarsQuadPt;
  
  PylithScalar density = 0;
  scalar_array properties(numPropsQuadPt);
  scalar_array stateVars(numVarsQuadPt);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(PylithScalar));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(PylithScalar));

    _matElastic->_calcDensity(&density, 
			      &properties[0], properties.size(),
			      &stateVars[0], stateVars.size());
    
    const PylithScalar densityE = data->density[iLoc];
    
    const PylithScalar tolerance = 1.0e-06;
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

  const int numLocs = data->numLocs;
  const int numPropsQuadPt = data->numPropsQuadPt;
  const int numVarsQuadPt = data->numVarsQuadPt;
  const int tensorSize = _matElastic->_tensorSize;
  
  scalar_array stress(tensorSize);
  scalar_array properties(numPropsQuadPt);
  scalar_array stateVars(numVarsQuadPt);
  scalar_array strain(tensorSize);
  scalar_array initialStress(tensorSize);
  scalar_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   properties.size()*sizeof(PylithScalar));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   stateVars.size()*sizeof(PylithScalar));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   strain.size()*sizeof(PylithScalar));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   initialStress.size()*sizeof(PylithScalar));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   initialStrain.size()*sizeof(PylithScalar));

    _matElastic->_calcStress(&stress[0], stress.size(),
			     &properties[0], properties.size(),
			     &stateVars[0], stateVars.size(),
			     &strain[0], strain.size(),
			     &initialStress[0], initialStress.size(),
			     &initialStrain[0], initialStrain.size(),
			     computeStateVars);

    const PylithScalar* stressE = &data->stress[iLoc*tensorSize];
    CPPUNIT_ASSERT(0 != stressE);

    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
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
      numConsts = 9;
      tensorSize = 3;
      break;
    case 3 :
      numConsts = 36;
      tensorSize = 6;
      break;
    } // switch
  const int numLocs = data->numLocs;
  const int numPropsQuadPt = data->numPropsQuadPt;
  const int numVarsQuadPt = data->numVarsQuadPt;
  
  scalar_array elasticConsts(numConsts);
  scalar_array properties(numPropsQuadPt);
  scalar_array stateVars(numVarsQuadPt);
  scalar_array strain(tensorSize);
  scalar_array initialStress(tensorSize);
  scalar_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(PylithScalar));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(PylithScalar));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));

    _matElastic->_calcElasticConsts(&elasticConsts[0], elasticConsts.size(),
				    &properties[0], properties.size(),
				    &stateVars[0], stateVars.size(),
				    &strain[0], strain.size(),
				    &initialStress[0], initialStress.size(),
				    &initialStrain[0], initialStrain.size());

    const PylithScalar* elasticConstsE = &data->elasticConsts[iLoc*numConsts];
    CPPUNIT_ASSERT(0 != elasticConstsE);
    
    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
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

  const int numLocs = data->numLocs;
  const int numPropsQuadPt = data->numPropsQuadPt;
  const int numVarsQuadPt = data->numVarsQuadPt;
  const int tensorSize = _matElastic->_tensorSize;
  
  scalar_array properties(numPropsQuadPt);
  scalar_array stateVars(numVarsQuadPt);
  scalar_array strain(tensorSize);
  scalar_array initialStress(tensorSize);
  scalar_array initialStrain(tensorSize);

  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &data->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(PylithScalar));
    memcpy(&stateVars[0], &data->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(PylithScalar));
    memcpy(&strain[0], &data->strain[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));
    memcpy(&initialStress[0], &data->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));
    memcpy(&initialStrain[0], &data->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(PylithScalar));

    _matElastic->_updateStateVars(&stateVars[0], stateVars.size(),
				  &properties[0], properties.size(),
				  &strain[0], strain.size(),
				  &initialStress[0], initialStress.size(),
				  &initialStrain[0], initialStrain.size());
    
    const PylithScalar* stateVarsE = 
      (numVarsQuadPt > 0) ? &data->stateVarsUpdated[iLoc*numVarsQuadPt] : 0;
    CPPUNIT_ASSERT( (0 < numVarsQuadPt && 0 != stateVarsE) ||
		    (0 == numVarsQuadPt && 0 == stateVarsE) );

    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
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

  const PylithScalar dt =
    _matElastic->_stableTimeStepImplicit(data->properties, data->numPropsQuadPt,
					 data->stateVars, data->numVarsQuadPt);

  const PylithScalar dtE = data->dtStableImplicit;

  const PylithScalar tolerance = 1.0e-06;
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

// ----------------------------------------------------------------------
// Setup mesh and material.
void
pylith::materials::TestElasticMaterial::_initialize(
					  topology::Mesh* mesh,
					  ElasticStrain1D* material,
					  const ElasticStrain1DData* data)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != material);
  CPPUNIT_ASSERT(0 != data);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/line3.mesh");
  iohandler.read(mesh);

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

  // Setup quadrature
  feassemble::Quadrature<topology::Mesh> quadrature;
  feassemble::GeometryLine1D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 1;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 1;
  const PylithScalar basis[] = { 0.455, -0.122, 0.667, -0.122, 0.455, 0.667 };
  const PylithScalar basisDeriv[] = { 
    -1.07735027e+00,
    -7.73502692e-02,
    1.15470054e+00,
    7.73502692e-02,
    1.07735027e+00,
    -1.15470054e+00,
  };
  const PylithScalar quadPtsRef[] = { -0.577350269, 0.577350269 };
  const PylithScalar quadWts[] = { 1.0, 1.0  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);


  // Get cells associated with material
  const int materialId = 24;
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);

  // Compute geometry for cells
  quadrature.initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  quadrature.computeGeometry(*mesh, cells);
#endif

  spatialdata::spatialdb::SimpleDB db;
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename("data/matinitialize.spatialdb");
  db.ioHandler(&dbIO);
  db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  spatialdata::spatialdb::SimpleDB dbStress;
  spatialdata::spatialdb::SimpleIOAscii dbIOStress;
  dbIOStress.filename("data/matstress.spatialdb");
  dbStress.ioHandler(&dbIOStress);
  dbStress.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  spatialdata::spatialdb::SimpleDB dbStrain;
  spatialdata::spatialdb::SimpleIOAscii dbIOStrain;
  dbIOStrain.filename("data/matstrain.spatialdb");
  dbStrain.ioHandler(&dbIOStrain);
  dbStrain.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  
  material->dbProperties(&db);
  material->id(materialId);
  material->label("my_material");
  material->normalizer(normalizer);
  material->dbInitialStress(&dbStress);
  material->dbInitialStrain(&dbStrain);
  
  material->initialize(*mesh, &quadrature);
} // _initialize


// End of file 
