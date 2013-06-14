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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestSubMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestSubMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestSubMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 2;
      const int ncorners = 3;
      const int cells[] = {
	0, 1, 3,
	0, 3, 2,
      };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
      const char* label = "bc";
      const int groupSize = 3;
      const int groupVertices[] = {
	1, 2, 3
      };
      const int submeshVertices[] = {
	3, 4, 5
      };
    } // _TestSubMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSubMesh::testConstructor(void)
{ // testConstructor
  SubMesh mesh;
  CPPUNIT_ASSERT(mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubMesh::testConstructorMesh(void)
{ // testConstructorMesh
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int nvertices = _TestSubMesh::groupSize;
  CPPUNIT_ASSERT_EQUAL(size_t(nvertices), vertices->size());
  int iV = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshVertices[iV++], *v_iter);
} // testConstructorMesh

// ----------------------------------------------------------------------
// Test sieveMesh().
void
pylith::topology::TestSubMesh::testSieveMesh(void)
{ // testSieveMesh
  Mesh mesh2D;
  _buildMesh(&mesh2D);

  SubMesh mesh(mesh2D, _TestSubMesh::label);
  
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());
} // testSieveMesh

// ----------------------------------------------------------------------
// Test createSubMesh().
void
pylith::topology::TestSubMesh::testCreateSubMesh(void)
{ // testCreateSubMesh
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  
  SubMesh mesh;
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT(!mesh._mesh.isNull());
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh.dimension());
  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const int nvertices = _TestSubMesh::groupSize;
  CPPUNIT_ASSERT_EQUAL(size_t(nvertices), vertices->size());
  int iV = 0;
  for (Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != vertices->end();
       ++v_iter)
    CPPUNIT_ASSERT_EQUAL(_TestSubMesh::submeshVertices[iV++], *v_iter);
} // testCreateSubMesh

// ----------------------------------------------------------------------
// Test coordsys().
void
pylith::topology::TestSubMesh::testCoordsys(void)
{ // testCoordsys
  Mesh mesh2D;
  _buildMesh(&mesh2D);

  SubMesh mesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim, mesh.coordsys()->spaceDim());
} // testCoordsys

// ----------------------------------------------------------------------
// Test debug().
void
pylith::topology::TestSubMesh::testDebug(void)
{ // testDebug
  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(false, mesh.debug());

  mesh.debug(true);
  CPPUNIT_ASSERT_EQUAL(true, mesh.debug());
} // testDebug

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::topology::TestSubMesh::testDimension(void)
{ // testDimension
  SubMesh mesh;
  CPPUNIT_ASSERT_EQUAL(0, mesh.dimension());

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh2(mesh2D, _TestSubMesh::label);
  CPPUNIT_ASSERT_EQUAL(_TestSubMesh::cellDim-1, mesh2.dimension());
} // testDimension

// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestSubMesh::testComm(void)
{ // testComm
  SubMesh mesh;

  Mesh mesh2D;
  _buildMesh(&mesh2D);
  mesh.createSubMesh(mesh2D, _TestSubMesh::label);

  CPPUNIT_ASSERT_EQUAL(PETSC_COMM_WORLD, mesh.comm());
} // testComm

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::topology::TestSubMesh::testInitialize(void)
{ // testInitialize
  Mesh mesh2D;
  _buildMesh(&mesh2D);
  SubMesh mesh(mesh2D, _TestSubMesh::label);

  mesh.initialize();
} // testInitialize

// ----------------------------------------------------------------------
void
pylith::topology::TestSubMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(0 != mesh);

  mesh->createSieveMesh(_TestSubMesh::cellDim);
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh->sieveMesh();

  ALE::Obj<Mesh::SieveMesh::sieve_type> sieve = 
    new Mesh::SieveMesh::sieve_type(sieveMesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  ALE::Obj<SieveFlexMesh::sieve_type> s = 
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  
  const int cellDim = _TestSubMesh::cellDim;
  const int ncells = _TestSubMesh::ncells;
  const int* cells = _TestSubMesh::cells;
  const int nvertices = _TestSubMesh::nvertices;
  const int ncorners = _TestSubMesh::ncorners;
  const int spaceDim = _TestSubMesh::cellDim;
  const PylithScalar* coordinates = _TestSubMesh::coordinates;
  const bool interpolate = false;
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, cellDim, ncells, (int*) cells,
					      nvertices, interpolate, 
					      ncorners);
  std::map<Mesh::SieveMesh::point_type,Mesh::SieveMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<Mesh::SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						       coordinates);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  typedef Mesh::IntSection::chart_type chart_type;
  const ALE::Obj<Mesh::IntSection>& groupField = 
    sieveMesh->getIntSection(_TestSubMesh::label);
  assert(!groupField.isNull());

  const int numPoints = _TestSubMesh::groupSize;
  const int numVertices = sieveMesh->depthStratum(0)->size();
  const int numCells = sieveMesh->heightStratum(0)->size();
  groupField->setChart(chart_type(numCells, numCells+numVertices));
  for(int i=0; i < numPoints; ++i)
    groupField->setFiberDimension(numCells+_TestSubMesh::groupVertices[i], 1);
  sieveMesh->allocate(groupField);
  
} // _buildMesh


// End of file 
