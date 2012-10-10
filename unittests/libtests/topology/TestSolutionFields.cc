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
  DM dmMesh = mesh.dmMesh();
  PetscErrorCode err;
  CPPUNIT_ASSERT(dmMesh);

  PetscInt       vStart, vEnd;
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  const char* labels[] = { "field A", "field B", "field C" };
  const int size = 3;
  const int fiberDimA = 2;
  const int fiberDimB = 3;
  const int fiberDimC = 4;

  for (int i=0; i < size; ++i)
    manager.add(labels[i], "displacement");

  Field<Mesh>& fieldA = manager.get(labels[0]);
  Field<Mesh>& fieldB = manager.get(labels[1]);
  Field<Mesh>& fieldC = manager.get(labels[2]);
  fieldA.newSection(FieldBase::VERTICES_FIELD, fiberDimA);

  fieldB.newSection(fieldA, fiberDimB);
  fieldC.newSection(fieldB, fiberDimC);

  manager.solutionName(labels[1]);
  const Field<Mesh>& solution = manager.solution();
  PetscSection section = solution.petscSection();
  CPPUNIT_ASSERT(section);
  PetscInt dof;
  err = PetscSectionGetDof(section, vStart, &dof);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(fiberDimB, dof);
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
