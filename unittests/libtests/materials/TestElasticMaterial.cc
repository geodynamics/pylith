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

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "data/ElasticMaterialData.hh" // USES ElasticMaterialData
#include "data/ElasticIsotropic3DData.hh" // USES ElasticIsotropic3DData

#include "pylith/feassemble/ParameterManager.hh" // USES ParameterManager

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
// Test density()
void
pylith::materials::TestElasticMaterial::testDensity(void)
{ // testDensity
  const double densityE[] = { 1.0, 8.0 };
  const int numQuadPts = 2;

  ElasticIsotropic3D material;
  material._density = (numQuadPts > 0) ? new double[numQuadPts] : 0;
  CPPUNIT_ASSERT(0 != material._density);
  memcpy(material._density, densityE, numQuadPts*sizeof(double));
  
  const double* density = material.density();
  for (int i=0; i < numQuadPts; ++i)
    CPPUNIT_ASSERT_EQUAL(densityE[i], density[i]);
} // testDensity

// ----------------------------------------------------------------------
// Test numElasticConsts() and elasticConsts()
void
pylith::materials::TestElasticMaterial::testElasticConsts(void)
{ // testElasticConsts
  const double elasticConstsE[] = { 4.0, 1.0, 2.0 };
  const int numQuadPts = 3;

  ElasticIsotropic3D material;
  material._elasticConsts = (numQuadPts > 0) ? new double[numQuadPts] : 0;
  CPPUNIT_ASSERT(0 != material._elasticConsts);
  memcpy(material._elasticConsts, elasticConstsE, numQuadPts*sizeof(double));
  
  const double* elasticConsts = material.elasticConsts();
  for (int i=0; i < numQuadPts; ++i)
    CPPUNIT_ASSERT_EQUAL(elasticConstsE[i], elasticConsts[i]);

  CPPUNIT_ASSERT_EQUAL(21, material.numElasticConsts());
} // testElasticConsts

// ----------------------------------------------------------------------
// Test calcProperties()
void
pylith::materials::TestElasticMaterial::testCalcProperties(void)
{ // testCalcProperties
  typedef ALE::Field::Mesh Mesh;
  typedef Mesh::sieve_type sieve_type;
  typedef Mesh::real_section_type real_section_type;

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
    ALE::New::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    mesh->setSieve(sieve);
    mesh->stratify();
    ALE::New::SieveBuilder<Mesh>::buildCoordinatesNew(mesh, spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  delete material._parameters; 
  material._parameters = new feassemble::ParameterManager(mesh);
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

  material.calcProperties(*cellIter, numQuadPts);

  const double tolerance = 1.0e-06;

  { // check density
  const double* density = material.density();
  const double* densityE = data.density;
  CPPUNIT_ASSERT(0 != density);
  CPPUNIT_ASSERT(0 != densityE);
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
  } // check density
  
  { // check elasticConsts
  const double* elasticConsts = material.elasticConsts();
  const double* elasticConstsE = data.elasticConsts;
  CPPUNIT_ASSERT(0 != elasticConsts);
  CPPUNIT_ASSERT(0 != elasticConstsE);
  const double size = numQuadPts * material.numElasticConsts();
  for (int i=0; i < size; ++i)
    if (fabs(elasticConstsE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i],
				   tolerance);
  } // check elasticConsts
} // testCalcProperties
    
// ----------------------------------------------------------------------
// Test _initCellData()
void
pylith::materials::TestElasticMaterial::testInitCellData(void)
{ // testInitCellData
  ElasticIsotropic3D material;

  CPPUNIT_ASSERT(0 == material._density);
  CPPUNIT_ASSERT(0 == material._elasticConsts);
  const int numQuadPts = 4;
  material._initCellData(numQuadPts);
  CPPUNIT_ASSERT(0 != material._density);
  CPPUNIT_ASSERT(0 != material._elasticConsts);
} // testInitCellData

// ----------------------------------------------------------------------
// Test _calcDensity()
void
pylith::materials::TestElasticMaterial::_testCalcDensity(
					ElasticMaterial* material,
					const ElasticMaterialData& data) const
{ // _testCalcDensity
  const int numQuadPts = data.numLocs;
  material->_initCellData(numQuadPts);
  material->_calcDensity(data.parameterData, data.numParameters, data.numLocs);
  const double* density = material->density();
  const double* densityE = data.density;
  
  CPPUNIT_ASSERT(0 != density);
  CPPUNIT_ASSERT(0 != densityE);

  const double tolerance = 1.0e-06;
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcElasticConsts()
void
pylith::materials::TestElasticMaterial::_testCalcElasticConsts(
				       ElasticMaterial* material,
				       const ElasticMaterialData& data) const
{ // _testCalcElasticConsts
  const int numQuadPts = data.numLocs;
  material->_initCellData(numQuadPts);
  material->_calcElasticConsts(data.parameterData, data.numParameters, 
			       data.numLocs);
  const double* elasticConsts = material->elasticConsts();
  const double* elasticConstsE = data.elasticConsts;
  
  CPPUNIT_ASSERT(0 != elasticConsts);
  CPPUNIT_ASSERT(0 != elasticConstsE);

  const double tolerance = 1.0e-06;
  const double size = numQuadPts * material->numElasticConsts();
  for (int i=0; i < size; ++i)
    if (fabs(elasticConstsE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, elasticConsts[i]/elasticConstsE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(elasticConstsE[i], elasticConsts[i],
				   tolerance);
} // _testCalcElasticConsts


// End of file 
