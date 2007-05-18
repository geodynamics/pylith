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
// Test clone()
void
pylith::materials::TestElasticMaterial::testClone(void)
{ // testClone
  const char* label = "my_material";
  const int id = 34;
  ElasticIsotropic3D material;
  material.label(label);
  material.id(id);

  ElasticMaterial* mycopy = material.clone();
  CPPUNIT_ASSERT(0 != mycopy);
  CPPUNIT_ASSERT(0 == strcmp(label, material.label().c_str()));
  CPPUNIT_ASSERT_EQUAL(id, material.id());

  delete mycopy; mycopy = 0;
} // testClone

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
    ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
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
  delete material._parameters; 
  material._parameters = new topology::FieldsManager(mesh);
  const int numQuadPts = 2;
  const int fiberDim = numQuadPts; // number of values in field per cell

  Mesh::label_sequence::iterator cellIter=cells->begin();

  material._parameters->addReal("density");
  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  parameterDensity->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterDensity);
  double cellData[numQuadPts];
  cellData[0] = data.parameterData[0];
  cellData[1] = data.parameterData[3];
  parameterDensity->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("mu");
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  parameterMu->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterMu);
  cellData[0] = data.parameterData[1];
  cellData[1] = data.parameterData[4];
  parameterMu->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("lambda");
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  parameterLambda->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterLambda);
  cellData[0] = data.parameterData[2];
  cellData[1] = data.parameterData[5];
  parameterLambda->updateAddPoint(*cellIter, cellData);

  material.initCellData(*cellIter, numQuadPts);
  const std::vector<double_array>& density = material.calcDensity();

  const double tolerance = 1.0e-06;
  const double* densityE = data.density;
  CPPUNIT_ASSERT(0 != densityE);
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i][0]/densityE[i], tolerance);
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
    ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
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
  delete material._parameters; 
  material._parameters = new topology::FieldsManager(mesh);
  const int numQuadPts = 2;
  const int fiberDim = numQuadPts; // number of values in field per cell

  Mesh::label_sequence::iterator cellIter=cells->begin();

  material._parameters->addReal("density");
  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  parameterDensity->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterDensity);
  double cellData[numQuadPts];
  cellData[0] = data.parameterData[0];
  cellData[1] = data.parameterData[3];
  parameterDensity->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("mu");
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  parameterMu->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterMu);
  cellData[0] = data.parameterData[1];
  cellData[1] = data.parameterData[4];
  parameterMu->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("lambda");
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  parameterLambda->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterLambda);
  cellData[0] = data.parameterData[2];
  cellData[1] = data.parameterData[5];
  parameterLambda->updateAddPoint(*cellIter, cellData);

  int strainSize = 0;
  switch(data.dimension)
    { // switch
    case 1 :
      strainSize = 1;
      break;
    case 2 :
      strainSize = 3;
      break;
    case 3 :
      strainSize = 6;
      break;
    } // switch
  std::vector<double_array> strain(numQuadPts);
  for (int iQuad=0, i=0; iQuad < numQuadPts; ++iQuad) {
    strain[iQuad].resize(strainSize);
    for (int iStrain=0; iStrain < strainSize; ++iStrain, ++i)
      strain[iQuad][iStrain] = data.strain[i];
  } // for

  material.initCellData(*cellIter, numQuadPts);
  const std::vector<double_array>& stress = material.calcStress(strain);

  const double tolerance = 1.0e-06;
  const double* stressE = data.stress;
  CPPUNIT_ASSERT(0 != stressE);
  
  for (int iQuad=0, i=0; iQuad < numQuadPts; ++iQuad) {
    const int stressSize = stress[iQuad].size();
    for (int iStress=0; iStress < stressSize; ++iStress, ++i) {
      if (fabs(stressE[i]) > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[iQuad][iStress]/stressE[i], 
				     tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[iQuad][iStress], 
				   tolerance);
    } // for
  } // for
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
    ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
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
  delete material._parameters; 
  material._parameters = new topology::FieldsManager(mesh);
  const int numQuadPts = 2;
  const int fiberDim = numQuadPts; // number of values in field per cell

  Mesh::label_sequence::iterator cellIter=cells->begin();

  material._parameters->addReal("density");
  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  parameterDensity->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterDensity);
  double cellData[numQuadPts];
  cellData[0] = data.parameterData[0];
  cellData[1] = data.parameterData[3];
  parameterDensity->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("mu");
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  parameterMu->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterMu);
  cellData[0] = data.parameterData[1];
  cellData[1] = data.parameterData[4];
  parameterMu->updateAddPoint(*cellIter, cellData);

  material._parameters->addReal("lambda");
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  parameterLambda->setFiberDimension(cells, fiberDim);
  mesh->allocate(parameterLambda);
  cellData[0] = data.parameterData[2];
  cellData[1] = data.parameterData[5];
  parameterLambda->updateAddPoint(*cellIter, cellData);

  int strainSize = 0;
  switch(data.dimension)
    { // switch
    case 1 :
      strainSize = 1;
      break;
    case 2 :
      strainSize = 3;
      break;
    case 3 :
      strainSize = 6;
      break;
    } // switch
  std::vector<double_array> strain(numQuadPts);
  for (int iQuad=0, i=0; iQuad < numQuadPts; ++iQuad) {
    strain[iQuad].resize(strainSize);
    for (int iStrain=0; iStrain < strainSize; ++iStrain, ++i)
      strain[iQuad][iStrain] = data.strain[i];
  } // for

  material.initCellData(*cellIter, numQuadPts);
  const std::vector<double_array>& elasticConsts = 
    material.calcDerivElastic(strain);

  const double tolerance = 1.0e-06;
  const double* elasticConstsE = data.elasticConsts;
  CPPUNIT_ASSERT(0 != elasticConstsE);

  for (int iQuad=0, i=0; iQuad < numQuadPts; ++iQuad) {
    const int numElasticConsts = elasticConsts[iQuad].size();
    for (int iConst=0; iConst < numElasticConsts; ++iConst, ++i) {
    if (fabs(elasticConstsE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
			   elasticConsts[iQuad][iConst]/elasticConstsE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], 
				   elasticConsts[iQuad][iConst],
				   tolerance);
    } // for
  } // for
} // testCalcDerivElastic
    
// ----------------------------------------------------------------------
// Test _calcDensity()
void
pylith::materials::TestElasticMaterial::_testCalcDensity(
					ElasticMaterial* material,
					const ElasticMaterialData& data) const
{ // _testCalcDensity
  const int numParameters = data.numParameters;
  std::vector<double_array> parameters(numParameters);
  for (int iParam=0, i=0; iParam < numParameters; ++iParam) {
    const int numValues = data.numParamValues[iParam];
    parameters[iParam].resize(numValues);
    for (int iValue=0; iValue < numValues; ++iValue)
      parameters[iParam][iValue] = data.parameterData[i++];
  } // for

  double_array density(1);
  material->_calcDensity(&density, parameters);
  const double* densityE = data.density;
  
  CPPUNIT_ASSERT(0 != densityE);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[0]/densityE[0], tolerance);
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcStress()
void
pylith::materials::TestElasticMaterial::_testCalcStress(
				       ElasticMaterial* material,
				       const ElasticMaterialData& data) const
{ // _testCalcElasticConsts
  const int numParameters = data.numParameters;
  std::vector<double_array> parameters(numParameters);
  for (int iParam=0, i=0; iParam < numParameters; ++iParam) {
    const int numValues = data.numParamValues[iParam];
    parameters[iParam].resize(numValues);
    for (int iValue=0; iValue < numValues; ++iValue)
      parameters[iParam][iValue] = data.parameterData[i++];
  } // for

  int stressSize = 0;
  switch(data.dimension)
    { // switch
    case 1 :
      stressSize = 1;
      break;
    case 2 :
      stressSize = 3;
      break;
    case 3 :
      stressSize = 6;
      break;
    } // switch
  double_array stress(stressSize);
  double_array strain(stressSize);
  for (int i=0; i < stressSize; ++i)
    strain[i] = data.strain[i];
  material->_calcStress(&stress, parameters, strain);
  const double* stressE = data.stress;
  
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
pylith::materials::TestElasticMaterial::_testCalcElasticConsts(
				       ElasticMaterial* material,
				       const ElasticMaterialData& data) const
{ // _testCalcElasticConsts
  const int numParameters = data.numParameters;
  std::vector<double_array> parameters(numParameters);
  for (int iParam=0, i=0; iParam < numParameters; ++iParam) {
    const int numValues = data.numParamValues[iParam];
    parameters[iParam].resize(numValues);
    for (int iValue=0; iValue < numValues; ++iValue)
      parameters[iParam][iValue] = data.parameterData[i++];
  } // for

  int numConsts = 0;
  int strainSize = 0;
  switch(data.dimension)
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
  double_array strain(strainSize);
  for (int i=0; i < strainSize; ++i)
    strain[i] = data.strain[i];
  material->_calcElasticConsts(&elasticConsts, parameters, strain);
  const double* elasticConstsE = data.elasticConsts;
  
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


// End of file 
