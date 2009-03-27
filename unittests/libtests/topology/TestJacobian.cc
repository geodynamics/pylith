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

#include "TestJacobian.hh" // Implementation of class methods

#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestJacobian );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestJacobian::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test matrix().
void
pylith::topology::TestJacobian::testMatrix(void)
{ // testMatrix
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);

  const PetscMat matrix = jacobian.matrix();
  CPPUNIT_ASSERT(0 != matrix);
} // testMatrix

// ----------------------------------------------------------------------
// Test assemble().
void
pylith::topology::TestJacobian::testAssemble(void)
{ // testAssemble
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);

  jacobian.assemble("flush_assembly");
  jacobian.assemble("final_assembly");
} // testAssemble

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestJacobian::testZero(void)
{ // testZero
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);

  jacobian.zero();
} // testZero

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestJacobian::testView(void)
{ // testView
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);

  jacobian.assemble("final_assembly");

  jacobian.view();
} // testView

// ----------------------------------------------------------------------
// Test write().
void
pylith::topology::TestJacobian::testWrite(void)
{ // testWrite
  Mesh mesh;
  SolutionFields fields(mesh);
  _initialize(&mesh, &fields);
  Jacobian jacobian(fields);

  jacobian.assemble("final_assembly");

  jacobian.write("jacobian.mat");
} // testWrite

// ----------------------------------------------------------------------
void
pylith::topology::TestJacobian::_initialize(Mesh* mesh,
					    SolutionFields* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  fields->add("disp t+dt", "displacement");
  fields->solutionName("disp t+dt");
  Field<Mesh>& solution = fields->solution();
  solution.newSection(FieldBase::VERTICES_FIELD, mesh->dimension());
  solution.allocate();
  solution.zero();
} // _initialize


// End of file 
