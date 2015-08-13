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

#include "TestElasticMaterial.hh" // Implementation of class methods

#include "data/ElasticMaterialData.hh" // USES ElasticMaterialData
#include "data/ElasticPlaneStrainData.hh" // USES ElasticPlaneStrainData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticMaterial );

// ----------------------------------------------------------------------
// Test dbInitialStress()
void
pylith::materials::TestElasticMaterial::testDBInitialStress(void)
{ // testDBInitialStress
  PYLITH_METHOD_BEGIN;

  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticPlaneStrain material;
  material.dbInitialStress(&db);
  
  CPPUNIT_ASSERT(material._dbInitialStress);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialStress->label()));

  PYLITH_METHOD_END;
} // testDBInitialStress

// ----------------------------------------------------------------------
// Test dbInitialStrain()
void
pylith::materials::TestElasticMaterial::testDBInitialStrain(void)
{ // testDBInitialStrain
  PYLITH_METHOD_BEGIN;

  const std::string& label = "my_database";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());
  
  ElasticPlaneStrain material;
  material.dbInitialStrain(&db);
  
  CPPUNIT_ASSERT(material._dbInitialStrain);
  CPPUNIT_ASSERT_EQUAL(label, std::string(material._dbInitialStrain->label()));

  PYLITH_METHOD_END;
} // testDBInitialStrain

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::materials::TestElasticMaterial::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  const PylithScalar tolerance = 1.0e-06;
  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Test initialStress field
  const PylithScalar *initialStressE = data.initialStress;
  CPPUNIT_ASSERT(initialStressE);
  int fiberDim = numQuadPts * tensorSize;

  CPPUNIT_ASSERT(material._initialFields);
  topology::Field& stressField = material._initialFields->get("initial stress");
  topology::VecVisitorMesh stressVisitor(stressField);
  const PetscScalar* initialStress = stressVisitor.localArray();
  PetscInt off = stressVisitor.sectionOffset(cell);
  CPPUNIT_ASSERT_EQUAL(fiberDim, stressVisitor.sectionDof(cell));
  for(int i = 0; i < fiberDim; ++i) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStress[off+i]/initialStressE[i]*data.pressureScale, tolerance);
  } // for

  // Test initialStrain field
  const PylithScalar *initialStrainE = data.initialStrain;
  CPPUNIT_ASSERT(initialStrainE);
  fiberDim = numQuadPts * tensorSize;

  CPPUNIT_ASSERT(material._initialFields);
  topology::Field& strainField = material._initialFields->get("initial strain");
  topology::VecVisitorMesh strainVisitor(strainField);
  const PetscScalar* initialStrain = strainVisitor.localArray();
  off = strainVisitor.sectionOffset(cell);
  CPPUNIT_ASSERT_EQUAL(fiberDim, strainVisitor.sectionDof(cell));
  for(int i = 0; i < fiberDim; ++i) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStrain[off+i]/initialStrainE[i], tolerance);
  } // for

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
      numElasticConsts = 9;
      break;
    case 3 :
      numElasticConsts = 21;
      break;
    default :
      assert(0);
    } // switch
  size = data.numLocs*numElasticConsts;
  CPPUNIT_ASSERT_EQUAL(size, material._elasticConstsCell.size());

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test retrievePropsAndVars().
void
pylith::materials::TestElasticMaterial::testRetrievePropsAndVars(void)
{ // testRetrievePropsAndVars
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();

  const PylithScalar tolerance = 1.0e-06;
  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;
  const int numVarsQuadPt = data.numVarsQuadPt;

  // Test cell arrays
  const PylithScalar* propertiesE = data.propertiesNondim;
  CPPUNIT_ASSERT(propertiesE);
  const scalar_array& properties = material._propertiesCell;
  size_t size = data.numLocs*data.numPropsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, properties.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, properties[i]/propertiesE[i],
				 tolerance);

  const PylithScalar* stateVarsE = data.stateVarsNondim;
  CPPUNIT_ASSERT( (0 < numVarsQuadPt && 0 != stateVarsE) ||
		  (0 == numVarsQuadPt && 0 == stateVarsE) );
  const scalar_array& stateVars = material._stateVarsCell;
  size = data.numLocs*numVarsQuadPt;
  CPPUNIT_ASSERT_EQUAL(size, stateVars.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i],
				 tolerance);

  const PylithScalar* initialStressE = data.initialStress;
  CPPUNIT_ASSERT(initialStressE);
  const scalar_array& initialStress = material._initialStressCell;
  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, initialStress.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStress[i]/initialStressE[i]*data.pressureScale,
				 tolerance);

  const PylithScalar* initialStrainE = data.initialStrain;
  CPPUNIT_ASSERT(initialStrainE);
  const scalar_array& initialStrain = material._initialStrainCell;
  size = data.numLocs*tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, initialStrain.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, initialStrain[i]/initialStrainE[i],
				 tolerance);

  PYLITH_METHOD_END;
} // testRetrievePropsAndVars

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticMaterial::testCalcDensity(void)
{ // testCalcDensity
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();
  const scalar_array& density = material.calcDensity();

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  const PylithScalar* densityE = data.density;
  CPPUNIT_ASSERT(densityE);
  const size_t size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, density.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i]*data.densityScale, tolerance);

  PYLITH_METHOD_END;
} // testCalcDensity
    
// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticMaterial::testCalcStress(void)
{ // testCalcStress
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Setup total strain
  scalar_array strain(data.strain, numQuadPts*tensorSize);

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();
  const scalar_array& stress = material.calcStress(strain);

  const PylithScalar* stressE = data.stress;
  CPPUNIT_ASSERT(stressE);
  const size_t size = numQuadPts * tensorSize;
  CPPUNIT_ASSERT_EQUAL(size, stress.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i]*data.pressureScale, tolerance);

  PYLITH_METHOD_END;
} // testCalcStress
    
// ----------------------------------------------------------------------
// Test calcDerivElastic()
void
pylith::materials::TestElasticMaterial::testCalcDerivElastic(void)
{ // testCalcDerivElastic
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  const int tensorSize = material._tensorSize;
  const int numQuadPts = data.numLocs;

  // Setup total strain
  scalar_array strain(data.strain, numQuadPts*tensorSize);

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();
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
  CPPUNIT_ASSERT(elasticConstsE);
  const size_t size = numQuadPts * numElasticConsts;
  CPPUNIT_ASSERT_EQUAL(size, elasticConsts.size());
  const PylithScalar tolerance = 1.0e-06;
  for (size_t i=0; i < size; ++i) {
    if (fabs(elasticConstsE[i]) > tolerance) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i]*data.pressureScale, tolerance);
    } else {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i]*data.pressureScale, tolerance);
    } // if/else
  } // for

  PYLITH_METHOD_END;
} // testCalcDerivElastic
    
// ----------------------------------------------------------------------
// Test updateStateVars()
void
pylith::materials::TestElasticMaterial::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  std::cout << "\n\nWARNING!! WARNING!! WARNING!!\n"
    "Need to implement using material with state variables.\n\n";

  PYLITH_METHOD_END;
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test calcStableTimeStepImplicit()
void
pylith::materials::TestElasticMaterial::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();
  const PylithScalar dt = material.stableTimeStepImplicit(mesh);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dtE = data.dtStableImplicit;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);

  PYLITH_METHOD_END;
} // testStableTimeStepImplicit

// ----------------------------------------------------------------------
// Test calcStableTimeStepExplicit()
void
pylith::materials::TestElasticMaterial::testStableTimeStepExplicit(void)
{ // testStableTimeStepExplicit
  PYLITH_METHOD_BEGIN;

  topology::Mesh mesh;
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _initialize(&mesh, &material, &data);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", materialId);
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
  PetscInt cell = cells[0];

  // Setup quadrature
  feassemble::Quadrature quadrature;
  feassemble::GeometryTri2D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 2;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const PylithScalar basis[numQuadPts*numCorners] = {
    1.0/6.0, 1.0/3.0, 1.0/2.0,
    1.0/6.0, 1.0/2.0, 1.0/3.0,
  };
  const PylithScalar basisDeriv[numQuadPts*numCorners*cellDim] = { 
    -0.5, 0.5,
    -0.5, 0.0,
     0.0, 0.5,
    -0.5, 0.5,
    -0.5, 0.0,
     0.0, 0.5,
  };
  const PylithScalar quadPtsRef[numQuadPts*spaceDim] = { 
    -1.0/3.0,      0.0,
         0.0, -1.0/3.0,
  };
  const PylithScalar quadWts[numQuadPts] = {
    1.0, 1.0,
  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  material.createPropsAndVarsVisitors();
  material.retrievePropsAndVars(cell);
  material.destroyPropsAndVarsVisitors();
  const PylithScalar dt = material.stableTimeStepExplicit(mesh, &quadrature);

  const PylithScalar tolerance = 1.0e-06;
  const PylithScalar dtE = 2.0*1.757359312880716 / 5196.15242;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);

  PYLITH_METHOD_END;
} // testStableTimeStepExplicit

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticMaterial::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestMaterial::setUp();
  _matElastic = 0;
  _dataElastic = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::materials::TestElasticMaterial::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  TestMaterial::tearDown();
  delete _matElastic; _matElastic = 0;
  delete _dataElastic; _dataElastic = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test _calcDensity()
void
pylith::materials::TestElasticMaterial::test_calcDensity(void)
{ // _testCalcDensity
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
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

  PYLITH_METHOD_END;
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcStress()
void
pylith::materials::TestElasticMaterial::test_calcStress(void)
{ // _testCalcStress
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
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
    CPPUNIT_ASSERT(stressE);

    const PylithScalar tolerance = (8 == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-04;
    for (int i=0; i < tensorSize; ++i)
      if (fabs(stressE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[i],
				     tolerance);
  } // for

  PYLITH_METHOD_END;
} // _testCalcStress

// ----------------------------------------------------------------------
// Test _calcElasticConsts()
void
pylith::materials::TestElasticMaterial::test_calcElasticConsts(void)
{ // _testCalcElasticConsts
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
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
    CPPUNIT_ASSERT(elasticConstsE);
    
    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for (int i=0; i < numConsts; ++i)
      if (fabs(elasticConstsE[i]) > tolerance) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i], 
				     tolerance);
      } else {
	const double stressScale = 1.0e+9;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i],
				     tolerance*stressScale);
      } // if/else
  } // for

  PYLITH_METHOD_END;
} // _testCalcElasticConsts

// ----------------------------------------------------------------------
// Test _updateStateVars()
void
pylith::materials::TestElasticMaterial::test_updateStateVars(void)
{ // test_updateStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
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

    const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
    for (int i=0; i < numVarsQuadPt; ++i)
      if (fabs(stateVarsE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], 
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i],
				     tolerance);
  } // for

  PYLITH_METHOD_END;
} // test_updateStateVars

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestElasticMaterial::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const PylithScalar dt =
    _matElastic->_stableTimeStepImplicit(data->properties, data->numPropsQuadPt,
					 data->stateVars, data->numVarsQuadPt);

  const PylithScalar dtE = data->dtStableImplicit;

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);

  PYLITH_METHOD_END;
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test _stableTimeStepExplicit()
void
pylith::materials::TestElasticMaterial::test_stableTimeStepExplicit(void)
{ // test_stableTimeStepExplicit
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_matElastic);
  CPPUNIT_ASSERT(_dataElastic);
  const ElasticMaterialData* data = _dataElastic;

  const PylithScalar minCellWidth = 1000.0;

  const PylithScalar dt =
    _matElastic->_stableTimeStepExplicit(data->properties, data->numPropsQuadPt,
					 data->stateVars, data->numVarsQuadPt,
					 minCellWidth);

  const PylithScalar dtE = data->dtStableExplicit;

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dtE, dt, tolerance); // TEMPORARY
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dt/dtE, tolerance);

  PYLITH_METHOD_END;
} // test_stableTimeStepExplicit

// ----------------------------------------------------------------------
// Setup nondimensionalization.
void
pylith::materials::TestElasticMaterial::setupNormalizer(void)
{ // setupNormalizer
  PYLITH_METHOD_BEGIN;

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.timeScale(_data->timeScale);
  normalizer.densityScale(_data->densityScale);
  _material->normalizer(normalizer);
  _matElastic->normalizer(normalizer);

  PYLITH_METHOD_END;
} // setupNormalizer

// ----------------------------------------------------------------------
// Setup mesh and material.
void
pylith::materials::TestElasticMaterial::_initialize(topology::Mesh* mesh,
						    ElasticPlaneStrain* material,
						    const ElasticPlaneStrainData* data)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(material);
  CPPUNIT_ASSERT(data);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  // Setup coordinates.
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Setup scales.
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(data->lengthScale);
  normalizer.pressureScale(data->pressureScale);
  normalizer.timeScale(data->timeScale);
  normalizer.densityScale(data->densityScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Setup quadrature
  feassemble::Quadrature quadrature;
  feassemble::GeometryTri2D geometry;
  quadrature.refGeometry(&geometry);
  const int cellDim = 2;
  const int numCorners = 3;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const PylithScalar basis[numQuadPts*numCorners] = {
    1.0/6.0, 1.0/3.0, 1.0/2.0,
    1.0/6.0, 1.0/2.0, 1.0/3.0,
  };
  const PylithScalar basisDeriv[numQuadPts*numCorners*cellDim] = { 
    -0.5, 0.5,
    -0.5, 0.0,
     0.0, 0.5,
    -0.5, 0.5,
    -0.5, 0.0,
     0.0, 0.5,
  };
  const PylithScalar quadPtsRef[numQuadPts*spaceDim] = { 
    -1.0/3.0,        0,
           0, -1.0/3.0,
  };
  const PylithScalar quadWts[numQuadPts] = {
    1.0, 1.0,
  };
  quadrature.initialize(basis, numQuadPts, numCorners,
			basisDeriv, numQuadPts, numCorners, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  // Get cells associated with material
  const int materialId = 24;
  PetscDM dmMesh = mesh->dmMesh();
  PetscIS cellIS = NULL;
  const PetscInt *cells = NULL;
  PetscInt numCells;
  PetscErrorCode err = 0;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", materialId, &cellIS);PYLITH_CHECK_ERROR(err);
  err = ISGetSize(cellIS, &numCells);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(cellIS, &cells);PYLITH_CHECK_ERROR(err);

  // Compute geometry for cells
  quadrature.initializeGeometry();
  err = ISRestoreIndices(cellIS, &cells);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&cellIS);PYLITH_CHECK_ERROR(err);

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

  PYLITH_METHOD_END;
} // _initialize


// End of file 
