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

#include "TestVertexFilterVecNorm.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/VertexFilterVecNorm.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Field.hh" // USES Field

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestVertexFilterVecNorm );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestVertexFilterVecNorm::testConstructor(void)
{ // testConstructor
  VertexFilterVecNorm<MeshField> filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestVertexFilterVecNorm::testFilter(void)
{ // testFilter
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef pylith::topology::Mesh::RealSection RealSection;

  const char* filename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int nvertices = 4;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::VECTOR;
  const PylithScalar fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
    3.1, 3.2,
    4.1, 4.2
  };
  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    topology::FieldBase::SCALAR;
  const int fiberDimE = 1;
  const PylithScalar fieldValuesE[] = {
    sqrt(pow(1.1, 2) + pow(1.2, 2)),
    sqrt(pow(2.1, 2) + pow(2.2, 2)),
    sqrt(pow(3.1, 2) + pow(3.2, 2)),
    sqrt(pow(4.1, 2) + pow(4.2, 2))
  };

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Set vertex field
  MeshField field(mesh);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());

  DM dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  PetscScalar *a;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  CPPUNIT_ASSERT_EQUAL(nvertices, vEnd-vStart);
  err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      a[off+d] = fieldValues[(v-vStart)*dof+d];
    }
  } // for
  err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);

  VertexFilterVecNorm<MeshField> filter;
  const MeshField& fieldF = filter.filter(field);
  PetscSection sectionF = fieldF.petscSection();
  Vec          vecF     = fieldF.localVector();
  CPPUNIT_ASSERT(sectionF);CPPUNIT_ASSERT(vecF);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(vecF, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(sectionF, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sectionF, v, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);
    for(PetscInt d = 0; d < dof; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a[off+d]/fieldValuesE[(v-vStart)*dof+d], tolerance);
    }
  } // for
  err = VecRestoreArray(vecF, &a);CHECK_PETSC_ERROR(err);
} // testFilter


// End of file 
