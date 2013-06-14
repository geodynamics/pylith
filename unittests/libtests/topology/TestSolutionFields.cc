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

#include "TestSolutionFields.hh" // Implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestSolutionFields );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSolutionFields::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test solutionName().
void
pylith::topology::TestSolutionFields::testSolutionName(void)
{ // testSolutionName
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const std::string& name = "my solution";
  manager.add(name.c_str(), "displacement");
  manager.solutionName(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, manager._solutionName);
} // testSolutionName

// ----------------------------------------------------------------------
// Test solution().
void
pylith::topology::TestSolutionFields::testSolution(void)
{ // testSolution
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const char* labels[] = { "field A", "field B", "field C" };
  const int size = 3;
  const int fiberDimA = 2;
  const int fiberDimB = 3;
  const int fiberDimC = 4;

  for (int i=0; i < size; ++i)
    manager.add(labels[i], "displacement");

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  Field<Mesh>& fieldA = manager.get(labels[0]);
  Field<Mesh>& fieldB = manager.get(labels[1]);
  Field<Mesh>& fieldC = manager.get(labels[2]);
  fieldA.newSection(vertices, fiberDimA);

  fieldB.newSection(fieldA, fiberDimB);
  fieldC.newSection(fieldB, fiberDimC);

  manager.solutionName(labels[1]);
  const Field<Mesh>& solution = manager.solution();
  const ALE::Obj<Mesh::RealSection>& sectionSoln = solution.section();
  CPPUNIT_ASSERT_EQUAL(fiberDimB,
		       sectionSoln->getFiberDimension(*(vertices->begin())));
} // testSolution

// ----------------------------------------------------------------------
void
pylith::topology::TestSolutionFields::_initialize(Mesh* mesh) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);
} // _initialize


// End of file 
