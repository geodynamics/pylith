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

#include "TestBoundaryMesh.hh" // Implementation of class methods

#include <Selection.hh> // USES submesh algorithms

#include "data/BoundaryMeshData.hh" // USES BoundaryMeshData

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMesh::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestBoundaryMesh::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test submesh() without fault.
void
pylith::bc::TestBoundaryMesh::testSubmesh(void)
{ // testSubmesh
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<Mesh> mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  const char* label = _data->bcLabel;
  const ALE::Obj<Mesh>& subMesh = 
    ALE::Selection<Mesh>::submesh(mesh, mesh->getIntSection(label));
  CPPUNIT_ASSERT(!subMesh.isNull());

  //subMesh->view("SUBMESH WITHOUT FAULT");

  const ALE::Obj<Mesh::label_sequence>& vertices = subMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, int(vertices->size()));

  int ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesNoFault[ipt], *v_iter);

  const ALE::Obj<Mesh::label_sequence>& cells = subMesh->heightStratum(1);
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = subMesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<Mesh> SieveAlg;

  CPPUNIT_ASSERT_EQUAL(_data->numCells, int(cells->size()));

  int icell = 0;
  int index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    const int depth = 1;
    const ALE::Obj<SieveAlg::coneArray>& cone = 
      SieveAlg::nCone(subMesh, *c_iter, depth);
    assert(!cone.isNull());
    const SieveAlg::coneArray::iterator vEnd = cone->end();

    CPPUNIT_ASSERT_EQUAL(_data->numCorners, int(cone->size()));

    int ivertex = 0;
    for(SieveAlg::coneArray::iterator v_iter=cone->begin();
	v_iter != vEnd;
	++v_iter, ++ivertex, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsNoFault[index], *v_iter);
  } // for
} // testSubmesh

// ----------------------------------------------------------------------
// Test submesh() with fault.
void
pylith::bc::TestBoundaryMesh::testSubmeshFault(void)
{ // testSubmeshFault
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<Mesh> mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  faults::FaultCohesiveKin fault;
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(mesh);

  const char* label = _data->bcLabel;
  const ALE::Obj<Mesh>& subMesh = 
    ALE::Selection<Mesh>::submesh(mesh, mesh->getIntSection(label));
  CPPUNIT_ASSERT(!subMesh.isNull());

  //subMesh->view("SUBMESH WITH FAULT");

  const ALE::Obj<Mesh::label_sequence>& vertices = subMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT_EQUAL(_data->numVerticesFault, int(vertices->size()));

  int ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[ipt], *v_iter);
    
  const ALE::Obj<Mesh::label_sequence>& cells = subMesh->depthStratum(1);
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = subMesh->getSieve();
  assert(!sieve.isNull());
  const int depth = 1;
  typedef ALE::SieveAlg<Mesh> SieveAlg;

  CPPUNIT_ASSERT_EQUAL(_data->numCells, int(cells->size()));

  int icell = 0;
  int index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    const ALE::Obj<SieveAlg::coneArray>& cone = 
      SieveAlg::nCone(subMesh, *c_iter, depth);
    assert(!cone.isNull());
    const SieveAlg::coneArray::iterator vEnd = cone->end();

    CPPUNIT_ASSERT_EQUAL(_data->numCorners, int(cone->size()));

    int ivertex = 0;
    for(SieveAlg::coneArray::iterator v_iter=cone->begin();
	v_iter != vEnd;
	++v_iter, ++ivertex, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsFault[index], *v_iter);
  } // for
} // testSubmeshFault


// End of file 
