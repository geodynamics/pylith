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

#include "TestOutputSolnSubset.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnSubset.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputSolnSubset );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputSolnSubset::testConstructor(void)
{ // testConstructor
  OutputSolnSubset output;
} // testConstructor

// ----------------------------------------------------------------------
// Test label()
void
pylith::meshio::TestOutputSolnSubset::testLabel(void)
{ // testLabel
  OutputSolnSubset output;

  const char* label = "boundary";

  output.label(label);
  CPPUNIT_ASSERT(0 == strcmp(label, output._label.c_str()));
} // testLabel

// ----------------------------------------------------------------------
// Test subdomainMesh()
void
pylith::meshio::TestOutputSolnSubset::testSubdomainMesh(void)
{ // testSubdomainMesh
  const char* filename = "data/quad4.mesh";
  const char* label = "bc3";
  const int nvertices = 3;
  const int verticesE[] = { 2, 4, 6 };
  const int ncells = 2;
  const int ncorners = 2;
  const int cellsE[] = { 2, 4, 4, 6 };

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename("data/quad4.mesh");
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());
  mesh->getFactory()->clear();

  OutputSolnSubset output;
  output.label(label);
  const ALE::Obj<Mesh> submesh = output.subdomainMesh(mesh);

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = submesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));
  int ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(verticesE[ipt], *v_iter);

  // Check cells
  const ALE::Obj<Mesh::label_sequence>& cells = submesh->heightStratum(1);
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = submesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<Mesh> SieveAlg;

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  int icell = 0;
  int index = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    const int depth = 1;
    const ALE::Obj<SieveAlg::coneArray>& cone = 
      SieveAlg::nCone(submesh, *c_iter, depth);
    assert(!cone.isNull());
    const SieveAlg::coneArray::iterator vEnd = cone->end();

    CPPUNIT_ASSERT_EQUAL(ncorners, int(cone->size()));
    
    for(SieveAlg::coneArray::iterator v_iter=cone->begin();
	v_iter != vEnd;
	++v_iter, ++index)
      CPPUNIT_ASSERT_EQUAL(cellsE[index], *v_iter);
  } // for
} // testSubdomainMesh


// End of file 
