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

#include "TestIntegratorInertia.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorInertia.hh" // USES IntegratorInertia
#include "data/IntegratorData.hh" // USES IntegratorData

// ----------------------------------------------------------------------
namespace pylith {
  namespace feassemble {
    class _TestIntegratorInertia;
  } // feassemble
} // pylith

/// Helper class for TestIntegrator
class pylith::feassemble::_TestIntegratorInertia {

public :
  /** Setup mesh.
   *
   * @param data Integrator data
   */
  static 
  ALE::Obj<ALE::Mesh>
  _setupMesh(const IntegratorData& data);
}; // _TestIntegrator

// ----------------------------------------------------------------------
// Test integrateLumped()
void
pylith::feassemble::TestIntegratorInertia::_testIntegrateLumped(
					   IntegratorInertia* integrator,
					   const IntegratorData& data) const
{ // _testIntegrateLumped
  typedef ALE::Mesh::real_section_type real_section_type;
  typedef ALE::Mesh::topology_type topology_type;

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
  typedef ALE::Mesh::topology_type topology_type;
  typedef topology_type::sieve_type sieve_type;

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


// End of file 
