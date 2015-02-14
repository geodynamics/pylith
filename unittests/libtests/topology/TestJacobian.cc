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
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  Jacobian jacobianB(field, "baij");
  Jacobian jacobianC(field, "baij", true);
  Jacobian jacobianD(field, "sbaij");
  Jacobian jacobianE(field, "sbaij", true);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test constructor with subdomain field
void
pylith::topology::TestJacobian::testConstructorSubDomain(void)
{ // testConstructorSubDomain
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);

  Mesh submesh(mesh, "bc");
  Field subfield(submesh);
  subfield.newSection(FieldBase::VERTICES_FIELD, submesh.dimension());
  subfield.allocate();
  subfield.zero();

  Jacobian jacobian(subfield);

  PYLITH_METHOD_END;
} // testConstructorSubDomain

// ----------------------------------------------------------------------
// Test matrix().
void
pylith::topology::TestJacobian::testMatrix(void)
{ // testMatrix
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  const PetscMat matrix = jacobian.matrix();
  CPPUNIT_ASSERT(matrix);

  PYLITH_METHOD_END;
} // testMatrix

// ----------------------------------------------------------------------
// Test assemble().
void
pylith::topology::TestJacobian::testAssemble(void)
{ // testAssemble
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("flush_assembly");
  jacobian.assemble("final_assembly");

  PYLITH_METHOD_END;
} // testAssemble

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestJacobian::testZero(void)
{ // testZero
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.zero();

  PYLITH_METHOD_END;
} // testZero

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestJacobian::testView(void)
{ // testView
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("final_assembly");

  jacobian.view();

  PYLITH_METHOD_END;
} // testView

// ----------------------------------------------------------------------
// Test write().
void
pylith::topology::TestJacobian::testWrite(void)
{ // testWrite
  PYLITH_METHOD_BEGIN;

  Mesh mesh;
  _initializeMesh(&mesh);
  Field field(mesh);
  _initializeField(&mesh, &field);
  Jacobian jacobian(field);

  jacobian.assemble("final_assembly");

  jacobian.write("jacobian.mat", mesh.comm());

  PYLITH_METHOD_END;
} // testWrite

// ----------------------------------------------------------------------
void
pylith::topology::TestJacobian::_initializeMesh(Mesh* mesh) const
{ // _initializeMesh
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);

  PYLITH_METHOD_END;
} // _initializeMesh

// ----------------------------------------------------------------------
void
pylith::topology::TestJacobian::_initializeField(Mesh* mesh,
                                                 Field* field) const
{ // _initializeField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(field);

  field->newSection(FieldBase::VERTICES_FIELD, mesh->dimension());
  field->allocate();
  field->zero();

  PYLITH_METHOD_END;
} // _initializeField

// End of file 
