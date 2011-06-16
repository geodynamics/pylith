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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestJacobian.hh" // Implementation of class methods

#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
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
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  Jacobian jacobianB(field, "baij");
  Jacobian jacobianC(field, "baij", true);
  Jacobian jacobianD(field, "sbaij");
  Jacobian jacobianE(field, "sbaij", true);
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor with subdomain field
void
pylith::topology::TestJacobian::testConstructorSubDomain(void)
{ // testConstructorSubDomain
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);

  SubMesh submesh(mesh, "bc");
  Field<SubMesh> subfield(submesh);
  subfield.newSection(FieldBase::VERTICES_FIELD, submesh.dimension());
  subfield.allocate();
  subfield.zero();

  Jacobian jacobian(subfield);
} // testConstructorSubDomain

// ----------------------------------------------------------------------
// Test matrix().
void
pylith::topology::TestJacobian::testMatrix(void)
{ // testMatrix
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  const PetscMat matrix = jacobian.matrix();
  CPPUNIT_ASSERT(0 != matrix);
} // testMatrix

// ----------------------------------------------------------------------
// Test assemble().
void
pylith::topology::TestJacobian::testAssemble(void)
{ // testAssemble
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("flush_assembly");
  jacobian.assemble("final_assembly");
} // testAssemble

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestJacobian::testZero(void)
{ // testZero
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.zero();
} // testZero

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestJacobian::testView(void)
{ // testView
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("final_assembly");

  jacobian.view();
} // testView

// ----------------------------------------------------------------------
// Test write().
void
pylith::topology::TestJacobian::testWrite(void)
{ // testWrite
  Mesh mesh;
  Field<Mesh> field(mesh);
  _initialize(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("final_assembly");

  jacobian.write("jacobian.mat", mesh.comm());
} // testWrite

// ----------------------------------------------------------------------
void
pylith::topology::TestJacobian::_initialize(Mesh* mesh,
                                            Field<Mesh>* field) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != field);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  field->newSection(FieldBase::VERTICES_FIELD, mesh->dimension());
  field->allocate();
  field->zero();
} // _initialize


// End of file 
