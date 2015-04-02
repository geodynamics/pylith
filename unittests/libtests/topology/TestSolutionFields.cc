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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestSolutionFields.hh" // Implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestSolutionFields );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSolutionFields::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  PYLITH_METHOD_END;
} // testConstructor
 
// ----------------------------------------------------------------------
// Test solutionName().
void
pylith::topology::TestSolutionFields::testSolutionName(void)
{ // testSolutionName
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const std::string& name = "my solution";
  manager.add(name.c_str(), "displacement");
  manager.solutionName(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, manager._solutionName);

  PYLITH_METHOD_END;
} // testSolutionName

// ----------------------------------------------------------------------
// Test solution().
void
pylith::topology::TestSolutionFields::testSolution(void)
{ // testSolution
  PYLITH_METHOD_BEGIN;

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

  Field& fieldA = manager.get(labels[0]);
  Field& fieldB = manager.get(labels[1]);
  Field& fieldC = manager.get(labels[2]);
  fieldA.newSection(FieldBase::VERTICES_FIELD, fiberDimA);
  fieldA.allocate();
  
  fieldB.newSection(fieldA, fiberDimB);
  fieldB.allocate();
  fieldC.newSection(fieldB, fiberDimC);
  fieldC.allocate();

  manager.solutionName(labels[1]);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum verticesStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();

  VecVisitorMesh solutionVisitor(manager.solution());
  CPPUNIT_ASSERT_EQUAL(fiberDimB, solutionVisitor.sectionDof(vStart));

  PYLITH_METHOD_END;
} // testSolution

// ----------------------------------------------------------------------
void
pylith::topology::TestSolutionFields::_initialize(Mesh* mesh) const
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  PYLITH_METHOD_END;
} // _initialize


// End of file 
