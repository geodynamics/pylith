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
// Test calcDensity()
void
pylith::materials::TestElasticMaterial::testCalcDensity(void)
{ // testCalcDensity
  typedef ALE::Mesh::topology_type topology_type;
  typedef topology_type::sieve_type sieve_type;
  typedef ALE::Mesh::real_section_type real_section_type;

  ALE::Obj<ALE::Mesh> mesh;
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

    mesh = new ALE::Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
    ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

    const bool interpolate = false;
    ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    sieve->stratify();
    topology->setPatch(0, sieve);
    topology->stratify();
    mesh->setTopology(topology);
    ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		   mesh->getRealSection("coordinates"), spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Mesh::int_section_type::patch_type patch = 0;
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  delete material._parameters; 
  material._parameters = new feassemble::ParameterManager(mesh);
  const int numQuadPts = 2;
  const int fiberDim = numQuadPts; // number of values in field per cell

  topology_type::label_sequence::iterator cellIter=cells->begin();

  material._parameters->addReal("density");
  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  parameterDensity->setFiberDimension(patch, cells, fiberDim);
  parameterDensity->allocate();
  double cellData[numQuadPts];
  cellData[0] = data.parameterData[0];
  cellData[1] = data.parameterData[3];
  parameterDensity->updateAdd(patch, *cellIter, cellData);

  material._parameters->addReal("mu");
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  parameterMu->setFiberDimension(patch, cells, fiberDim);
  parameterMu->allocate();
  cellData[0] = data.parameterData[1];
  cellData[1] = data.parameterData[4];
  parameterMu->updateAdd(patch, *cellIter, cellData);

  material._parameters->addReal("lambda");
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  parameterLambda->setFiberDimension(patch, cells, fiberDim);
  parameterLambda->allocate();
  cellData[0] = data.parameterData[2];
  cellData[1] = data.parameterData[5];
  parameterLambda->updateAdd(patch, *cellIter, cellData);

  material._initCellData(numQuadPts);
  const double* density = material.calcDensity(*cellIter, patch, numQuadPts);

  const double tolerance = 1.0e-06;
  const double* densityE = data.density;
  CPPUNIT_ASSERT(0 != density);
  CPPUNIT_ASSERT(0 != densityE);
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
} // testCalcProperties
    
// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticMaterial::testCalcStress(void)
{ // testCalcProperties
  typedef ALE::Mesh::topology_type topology_type;
  typedef topology_type::sieve_type sieve_type;
  typedef ALE::Mesh::real_section_type real_section_type;

  ALE::Obj<ALE::Mesh> mesh;
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

    mesh = new ALE::Mesh(PETSC_COMM_WORLD, cellDim);
    ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
    ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

    const bool interpolate = false;
    ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, cellDim, numCells,
	       const_cast<int*>(cells), numVertices, interpolate, numCorners);
    sieve->stratify();
    topology->setPatch(0, sieve);
    topology->stratify();
    mesh->setTopology(topology);
    ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		   mesh->getRealSection("coordinates"), spaceDim, vertCoords);
  } // create mesh

  // Get cells associated with material
  const ALE::Mesh::int_section_type::patch_type patch = 0;
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<topology_type>& topology = coordinates->getTopology();
  const ALE::Obj<topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);

  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  delete material._parameters; 
  material._parameters = new feassemble::ParameterManager(mesh);
  const int numQuadPts = 2;
  const int fiberDim = numQuadPts; // number of values in field per cell

  topology_type::label_sequence::iterator cellIter=cells->begin();

  material._parameters->addReal("density");
  const ALE::Obj<real_section_type>& parameterDensity = 
    material._parameters->getReal("density");
  parameterDensity->setFiberDimension(patch, cells, fiberDim);
  parameterDensity->allocate();
  double cellData[numQuadPts];
  cellData[0] = data.parameterData[0];
  cellData[1] = data.parameterData[3];
  parameterDensity->updateAdd(patch, *cellIter, cellData);

  material._parameters->addReal("mu");
  const ALE::Obj<real_section_type>& parameterMu = 
    material._parameters->getReal("mu");
  parameterMu->setFiberDimension(patch, cells, fiberDim);
  parameterMu->allocate();
  cellData[0] = data.parameterData[1];
  cellData[1] = data.parameterData[4];
  parameterMu->updateAdd(patch, *cellIter, cellData);

  material._parameters->addReal("lambda");
  const ALE::Obj<real_section_type>& parameterLambda = 
    material._parameters->getReal("lambda");
  parameterLambda->setFiberDimension(patch, cells, fiberDim);
  parameterLambda->allocate();
  cellData[0] = data.parameterData[2];
  cellData[1] = data.parameterData[5];
  parameterLambda->updateAdd(patch, *cellIter, cellData);

  material._initCellData(numQuadPts);
  const double* stress = material.calcStress(*cellIter, patch, data.strain, 
					     numQuadPts, data.dimension);

  const double tolerance = 1.0e-06;
  const double* stressE = data.stress;
  CPPUNIT_ASSERT(0 != stress);
  CPPUNIT_ASSERT(0 != stressE);
  const double size = numQuadPts * material.stressSize();
  for (int i=0; i < size; ++i)
    if (fabs(stressE[i]) > tolerance)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stress[i]/stressE[i], 
				   tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stressE[i], stress[i], tolerance);
} // testCalcStress
    
// ----------------------------------------------------------------------
// Test calcDerivElastic()
void
pylith::materials::TestElasticMaterial::testCalcDerivElastic(void)
{ // testCalcDerivElastic
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

  material._initCellData(numQuadPts);
  const double* elasticConsts = 
    material.calcDerivElastic(*cellIter, patch, data.strain, 
			      numQuadPts, data.dimension);

  const double tolerance = 1.0e-06;
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
} // testCalcDerivElastic
    
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
  CPPUNIT_ASSERT(0 != material._stress);
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
  const double* densityE = data.density;
  const double* density = material->_density;
  
  CPPUNIT_ASSERT(0 != density);
  CPPUNIT_ASSERT(0 != densityE);

  const double tolerance = 1.0e-06;
  const double size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density[i]/densityE[i], tolerance);
} // _testCalcDensity

// ----------------------------------------------------------------------
// Test _calcStress()
void
pylith::materials::TestElasticMaterial::_testCalcStress(
				       ElasticMaterial* material,
				       const ElasticMaterialData& data) const
{ // _testCalcElasticConsts
  const int numQuadPts = data.numLocs;
  material->_initCellData(numQuadPts);
  material->_calcStress(data.parameterData, data.numParameters, data.strain, 
			data.numLocs, data.dimension);
  const double* stressE = data.stress;
  const double* stress = material->_stress;
  
  CPPUNIT_ASSERT(0 != stress);
  CPPUNIT_ASSERT(0 != stressE);

  const double tolerance = 1.0e-06;
  const double size = numQuadPts * material->stressSize();
  for (int i=0; i < size; ++i)
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
  const int numQuadPts = data.numLocs;
  material->_initCellData(numQuadPts);
  material->_calcElasticConsts(data.parameterData, data.numParameters, 
			       data.strain, data.numLocs, data.dimension);
  const double* elasticConstsE = data.elasticConsts;
  const double* elasticConsts = material->_elasticConsts;
  
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
