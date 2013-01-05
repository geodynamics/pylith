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

#include "TestCellFilterAvg.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/CellFilterAvg.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestCellFilterAvg );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestCellFilterAvg::testConstructor(void)
{ // testConstructor
  CellFilterAvg<topology::Mesh, MeshField> filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestCellFilterAvg::testFilter(void)
{ // testFilter
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;
  typedef pylith::topology::Mesh::RealSection RealSection;

  const char* filename = "data/quad4.mesh";
  const int fiberDim = 2;
  const int ncells = 2;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::MULTI_SCALAR;
  const PylithScalar fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
  };
  const int cellDim = 2;
  const int numBasis = 4;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const PylithScalar basis[] = {
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
  };
  const PylithScalar basisDerivRef[] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
  };
  const PylithScalar quadPtsRef[] = {
    1.0, 0.0,
   -1.0, 0.0,};
  const PylithScalar quadWts[] = { 1.5, 0.5 };
  const PylithScalar minJacobian = 1.0;

  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    topology::FieldBase::SCALAR;
  const int fiberDimE = 1;
  const PylithScalar fieldValuesE[] = {
    (1.5*1.1 + 0.5*1.2)/2.0,
    (1.5*2.1 + 0.5*2.2)/2.0,
  };

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Set cell field
  MeshField field(mesh);
  field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());

  DM dmMesh = mesh.dmMesh();
  PetscInt       cStart, cEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  PetscScalar *a;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  CPPUNIT_ASSERT_EQUAL(ncells, cEnd-cStart);
  err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, c, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, c, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      a[off+d] = fieldValues[(c-cStart)*dof+d];
    }
  } // for
  err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);

  feassemble::Quadrature<topology::Mesh> quadrature;
  quadrature.initialize(basis, numQuadPts, numBasis,
			basisDerivRef, numQuadPts, numBasis, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  CellFilterAvg<topology::Mesh, MeshField> filter;
  filter.quadrature(&quadrature);

  const topology::Field<topology::Mesh>& fieldF = filter.filter(field);
  PetscSection sectionF = fieldF.petscSection();
  Vec          vecF     = fieldF.localVector();
  CPPUNIT_ASSERT(sectionF);CPPUNIT_ASSERT(vecF);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));

  const PylithScalar tolerance = 1.0e-06;
  err = VecGetArray(vecF, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt dof, off;

    err = PetscSectionGetDof(sectionF, c, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(sectionF, c, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);
    for(PetscInt d = 0; d < dof; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, a[off+d]/fieldValuesE[(c-cStart)*dof+d], tolerance);
    }
  } // for
  err = VecRestoreArray(vecF, &a);CHECK_PETSC_ERROR(err);
} // testFilter


// End of file 
