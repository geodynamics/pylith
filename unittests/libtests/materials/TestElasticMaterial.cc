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
  const int numQuadPts = 2;
  const int_array& numParamValues = material._getNumParamValues();
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;
  
  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._parameters = new real_section_type(mesh->comm(), mesh->debug());
  material._parameters->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._parameters);

  material._parameters->updatePoint(*c_iter, data.parameterData);

  material.getStateVarsCell(*c_iter, numQuadPts);
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
  const int numQuadPts = 2;
  const int_array& numParamValues = material._getNumParamValues();
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._parameters = new real_section_type(mesh->comm(), mesh->debug());
  material._parameters->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._parameters);

  material._parameters->updatePoint(*c_iter, data.parameterData);

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
  double_array strain(data.strain, numQuadPts*strainSize);

  material.getStateVarsCell(*c_iter, numQuadPts);
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
  const int numQuadPts = 2;
  const int_array& numParamValues = material._getNumParamValues();
  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  Mesh::label_sequence::iterator c_iter = cells->begin();

  const int fiberDim = numQuadPts * numParamsQuadPt;
  material._parameters = new real_section_type(mesh->comm(), mesh->debug());
  material._parameters->setFiberDimension(cells, fiberDim);
  mesh->allocate(material._parameters);

  material._parameters->updatePoint(*c_iter, data.parameterData);

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
  double_array strain(data.strain, numQuadPts*strainSize);

  material.getStateVarsCell(*c_iter, numQuadPts);
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
// Test _calcDensity()
void
pylith::materials::TestElasticMaterial::_testCalcDensity(
					ElasticMaterial* material,
					const ElasticMaterialData& data) const
{ // _testCalcDensity
  double_array density(1);
  material->_calcDensity(&density[0], data.parameterData, data.numParamsQuadPt);

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
  material->_calcStress(&stress[0], stress.size(),
			data.parameterData, data.numParamsQuadPt,
			data.strain, stressSize);

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
  material->_calcElasticConsts(&elasticConsts[0], numConsts,
			       data.parameterData, data.numParamsQuadPt,
			       data.strain, strainSize);

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
