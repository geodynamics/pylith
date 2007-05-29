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

#include "TestElasticityExplicit.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityExplicit.hh" // USES ElasticityExplicit
//#include "data/ElasticityExplicitData.hh" // USES ElasticityExplicitData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestElasticityExplicit::testConstructor(void)
{ // testConstructor
  ElasticityExplicit integrator;
} // testConstructor

// ----------------------------------------------------------------------
// Test timeStep().
void
pylith::feassemble::TestElasticityExplicit::testTimeStep(void)
{ // testTimeStep
  ElasticityExplicit integrator;

  const double dt1 = 2.0;
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dtm1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);
} // testTimeStep

// ----------------------------------------------------------------------
// Test StableTimeStep().
void
pylith::feassemble::TestElasticityExplicit::testStableTimeStep(void)
{ // testStableTimeStep
  ElasticityExplicit integrator;

  const double dt1 = 2.0;
  integrator.timeStep(dt1);
  const double stableTimeStep = integrator.stableTimeStep();
  CPPUNIT_ASSERT_EQUAL(dt1, stableTimeStep);
} // testStableTimeStep

// ----------------------------------------------------------------------
// Test material().
void
pylith::feassemble::TestElasticityExplicit::testMaterial(void)
{ // testMaterial
  ElasticityExplicit integrator;

  materials::ElasticIsotropic3D material;
  const int id = 3;
  const std::string label("my material");
  material.id(id);
  material.label(label.c_str());
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(id, integrator._material->id());
  CPPUNIT_ASSERT_EQUAL(label, std::string(integrator._material->label()));
} // testMaterial

// ----------------------------------------------------------------------
// Test updateState().
void 
pylith::feassemble::TestElasticityExplicit::testUpdateState(void)
{ // testUpdateState
  throw std::logic_error("Unit test not implemented.");
} // testUpdateState

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateResidual(void)
{ // testIntegrateResidual
  throw std::logic_error("Unit test not implemented.");
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  throw std::logic_error("Unit test not implemented.");
} // testIntegrateJacobian


#if 0
// ----------------------------------------------------------------------
// Test integrateLumped()
void
pylith::feassemble::TestIntegratorInertia::_testIntegrateLumped(
					   IntegratorInertia* integrator,
					   const IntegratorData& data) const
{ // _testIntegrateLumped

  CPPUNIT_ASSERT(false);

  ALE::Obj<ALE::Mesh> mesh = _TestIntegratorInertia::_setupMesh(data);
  const ALE::Mesh::int_section_type::patch_type patch = 0;

  // Fiber dimension (number of values in field per vertex) for fields
  const int fiberDim = data.fiberDim;

  // Setup field for action result
  const ALE::Obj<real_section_type>& fieldOut =
    mesh->getRealSection("fieldOut");
  fieldOut->setName("fieldOut");
  fieldOut->setFiberDimensionByDepth(patch, 0, fiberDim);
  fieldOut->allocate();

  // Should read density from spatial database
  const ALE::Obj<real_section_type>& density =
    mesh->getRealSection("density");
  density->setName("density");
  density->setFiberDimensionByDepth(patch, 0, 
				    integrator->_quadrature->numQuadPts());
  density->allocate();
  integrator->setDensity(density);

  // Integrate action
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  integrator->integrateLumped(fieldOut, coordinates);
  //fieldOut->view("field out");
  
  // Check values in output field
  int iVertex = 0;
  const ALE::Obj<topology_type::label_sequence>& vertices = 
    mesh->getTopology()->depthStratum(patch, 0);
  const topology_type::label_sequence::iterator verticesEnd =
    vertices->end();
  const double tolerance = 1.0e-06;
  for (topology_type::label_sequence::iterator vIter=vertices->begin();
       vIter != verticesEnd;
       ++vIter, ++iVertex) {
    const real_section_type::value_type* vals = 
      fieldOut->restrict(patch, *vIter);
    const double* valsE = &data.valsLumped[iVertex*fiberDim];
    const int dim = fieldOut->getFiberDimension(patch, *vIter);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dim);
    for (int iDim=0; iDim < fiberDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[iDim]/valsE[iDim], tolerance);
  } // for
} // _testIntegrateLumped


// ----------------------------------------------------------------------
// Setup mesh.
ALE::Obj<ALE::Mesh>
pylith::feassemble::_TestIntegratorInertia::_setupMesh(const IntegratorData& data)
{ // _setupMesh
  const int cellDim = data.cellDim;
  const int numCorners = data.numCorners;
  const int spaceDim = data.spaceDim;
  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  CPPUNIT_ASSERT(0 != vertCoords);
  CPPUNIT_ASSERT(0 != cells);

  ALE::Obj<ALE::Mesh> mesh = new ALE::Mesh(PETSC_COMM_WORLD, cellDim);
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

  return mesh;
} // _setupMesh
#endif


// End of file 
