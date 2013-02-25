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

#include "TestBoundaryMesh.hh" // Implementation of class methods

#include "data/BoundaryMeshData.hh" // USES BoundaryMeshData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

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

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  mesh.nondimensionalize(normalizer);

  // Create submesh
  topology::SubMesh submesh(mesh, _data->bcLabel);
  CPPUNIT_ASSERT(submesh.dmMesh());

  // Check vertices
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    submesh.sieveMesh()->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, int(vertices->size()));

  int ipt = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesNoFault[ipt], *v_iter);

  // Check cells
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh.sieveMesh()->heightStratum(1);
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = 
    submesh.sieveMesh()->getSieve();
  assert(!sieve.isNull());
  CPPUNIT_ASSERT_EQUAL(_data->numCells, int(cells->size()));

  ALE::ISieveVisitor::NConeRetriever<SieveSubMesh::sieve_type> ncV(*sieve, (int) pow(sieve->getMaxConeSize(), submesh.sieveMesh()->depth()));

  int icell = 0;
  int index = 0;
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    ALE::ISieveTraversal<SieveSubMesh::sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    const SieveSubMesh::point_type *cone = ncV.getPoints();
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

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  spatialdata::units::Nondimensional normalizer;
  cs.setSpaceDim(mesh.dimension());
  cs.initialize();
  mesh.coordsys(&cs);
  mesh.nondimensionalize(normalizer);

  // Adjust topology
  faults::FaultCohesiveKin fault;
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = mesh.sieveMesh()->getIntSection(_data->faultLabel)->size();
  int firstFaultCell      = mesh.sieveMesh()->getIntSection(_data->faultLabel)->size();
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(&mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);

  // Create submesh
  topology::SubMesh submesh(mesh, _data->bcLabel);
  CPPUNIT_ASSERT(!submesh.sieveMesh().isNull());

  // Check vertices
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    submesh.sieveMesh()->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  CPPUNIT_ASSERT_EQUAL(_data->numVerticesFault, int(vertices->size()));

  int ipt = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt)
    CPPUNIT_ASSERT_EQUAL(_data->verticesFault[ipt], *v_iter);

  // Check cells
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh.sieveMesh()->heightStratum(1);
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = 
    submesh.sieveMesh()->getSieve();
  assert(!sieve.isNull());
  CPPUNIT_ASSERT_EQUAL(_data->numCells, int(cells->size()));

  ALE::ISieveVisitor::NConeRetriever<SieveSubMesh::sieve_type> ncV(*sieve, (int) pow(sieve->getMaxConeSize(), submesh.sieveMesh()->depth()));

  int icell = 0;
  int index = 0;
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++icell) {
    ALE::ISieveTraversal<SieveSubMesh::sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int coneSize = ncV.getSize();
    const SieveSubMesh::point_type *cone = ncV.getPoints();
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, coneSize);

    for(int v = 0; v < coneSize; ++v, ++index)
      CPPUNIT_ASSERT_EQUAL(_data->cellsFault[index], cone[v]);
    ncV.clear();
  } // for
} // testSubmeshFault


// End of file 
