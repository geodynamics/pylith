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

#include "data/BoundaryMeshData.hh" // USES BoundaryMeshData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include <Selection.hh> // USES submesh algorithms

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

  topology::Mesh mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const char* label = _data->bcLabel;
  const ALE::Obj<SieveSubMesh>& subMesh = 
    ALE::Selection<Mesh>::submeshV<SieveSubMesh>(sieveMesh, 
					    sieveMesh->getIntSection(label));
  CPPUNIT_ASSERT(!subMesh.isNull());

  //subMesh->view("SUBMESH WITHOUT FAULT");

  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = subMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, int(vertices->size()));

  int ipt = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesNoFault[ipt], *v_iter);

  const ALE::Obj<SieveSubMesh::label_sequence>& cells = subMesh->heightStratum(1);
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = subMesh->getSieve();
  assert(!sieve.isNull());

  CPPUNIT_ASSERT_EQUAL(_data->numCells, (int) cells->size());

  ALE::ISieveVisitor::NConeRetriever<sieve_type> ncV(*sieve, (int) pow(sieve->getMaxConeSize(), subMesh->depth()));

  int icell = 0;
  int index = 0;
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    ALE::ISieveTraversal<sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    const SieveMesh::point_type *cone = ncV.getPoints();

    CPPUNIT_ASSERT_EQUAL(_data->numCorners, coneSize);

    for(int v = 0; v < coneSize; ++v, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsNoFault[index], cone[v]);
    ncV.clear();
  } // for
} // testSubmesh

// ----------------------------------------------------------------------
// Test submesh() with fault.
void
pylith::bc::TestBoundaryMesh::testSubmeshFault(void)
{ // testSubmeshFault
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filename);
  iohandler.read(&mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  faults::FaultCohesiveKin fault;
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(sieveMesh, _flipFault);

  const char* label = _data->bcLabel;
  const ALE::Obj<SieveSubMesh>& subMesh = 
    ALE::Selection<SieveMesh>::submeshV<SieveSubMesh>(sievMesh,
						      sievMesh->getIntSection(label));
  CPPUNIT_ASSERT(!subMesh.isNull());

  //subMesh->view("Submesh for mesh w/fault");
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = subMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  CPPUNIT_ASSERT_EQUAL(_data->numVerticesFault, int(vertices->size()));

  int ipt = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[ipt], *v_iter);
    
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = subMesh->depthStratum(1);
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = subMesh->getSieve();
  assert(!sieve.isNull());
  const int depth = 1;
  typedef ALE::SieveAlg<SieveMesh> SieveAlg;

  CPPUNIT_ASSERT_EQUAL(_data->numCells, int(cells->size()));

  ALE::ISieveVisitor::NConeRetriever<sieve_type> ncV(*sieve, (int) pow(sieve->getMaxConeSize(), subMesh->depth()));

  int icell = 0;
  int index = 0;
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    ALE::ISieveTraversal<sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    const SieveMesh::point_type *cone = ncV.getPoints();

    CPPUNIT_ASSERT_EQUAL(_data->numCorners, coneSize);

    for(int v = 0; v < coneSize; ++v, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsFault[index], cone[v]);
    ncV.clear();
  } // for
} // testSubmeshFault


// End of file 
