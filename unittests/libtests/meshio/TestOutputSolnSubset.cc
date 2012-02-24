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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestOutputSolnSubset.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnSubset.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputSolnSubset );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

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

  topology::Mesh mesh;
  MeshIOAscii iohandler;
  iohandler.filename("data/quad4.mesh");
  iohandler.read(&mesh);

  OutputSolnSubset output;
  output.label(label);
  const topology::SubMesh& submesh = output.subdomainMesh(mesh);
  const ALE::Obj<SieveMesh>& sieveSubMesh = submesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveSubMesh.isNull());

  // Check vertices
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveSubMesh->depthStratum(0);
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));
  int ipt = 0;
  for (SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(verticesE[ipt], *v_iter);

  // Check cells
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = sieveSubMesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<SieveMesh> SieveAlg;

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  ALE::ISieveVisitor::PointRetriever<SieveSubMesh::sieve_type> pV(sieve->getMaxConeSize());
  int i = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin(); c_iter != cellsEnd; ++c_iter) {
    sieve->cone(*c_iter, pV);
    const SieveSubMesh::point_type *cone = pV.getPoints();
    CPPUNIT_ASSERT_EQUAL(ncorners, (int) pV.getSize());
    for(int p = 0; p < pV.getSize(); ++p, ++i) {
      CPPUNIT_ASSERT_EQUAL(cellsE[i], cone[p]);
    }
    pV.clear();
  } // for
} // testSubdomainMesh


// End of file 
