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

#include "TestCellFilterAvg.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/CellFilterAvg.hh"

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestCellFilterAvg );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestCellFilterAvg::testConstructor(void)
{ // testConstructor
  CellFilterAvg filter;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter() w/mesh.
void
pylith::meshio::TestCellFilterAvg::testFilterMesh(void)
{ // testFilterMesh
  PYLITH_METHOD_BEGIN;

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
  const PylithScalar fieldScale = 2.0;

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  // Set cell field
  topology::Field field(mesh);
  field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());
  field.scale(fieldScale);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    const PetscInt off = fieldVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(c));

    for(PetscInt d = 0; d < fiberDim; ++d, ++index) {
      fieldArray[off+d] = fieldValues[index];
    } // for
  } // for

  feassemble::Quadrature quadrature;
  quadrature.initialize(basis, numQuadPts, numBasis,
			basisDerivRef, numQuadPts, numBasis, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  CellFilterAvg filter;
  filter.quadrature(&quadrature);

  const topology::Field& fieldF = filter.filter(field);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));
  CPPUNIT_ASSERT_EQUAL(fieldScale, fieldF.scale());

  topology::VecVisitorMesh fieldFVisitor(fieldF);
  const PetscScalar* fieldFArray = fieldFVisitor.localArray();CPPUNIT_ASSERT(fieldFArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    const PetscInt off = fieldFVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldFVisitor.sectionDof(c));
    for(PetscInt d = 0; d < fiberDimE; ++d, ++index) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldFArray[off+d]/fieldValuesE[index], tolerance);
    }
  } // for

  PYLITH_METHOD_END;
} // testFilterMesh


// ----------------------------------------------------------------------
// Test filter() w/submesh.
void
pylith::meshio::TestCellFilterAvg::testFilterSubMesh(void)
{ // testFilterMesh
  PYLITH_METHOD_BEGIN;

  const char* filename = "data/quad4.mesh";
  const char* group = "bc";

  const int fiberDim = 2;
  const int ncells = 2;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::MULTI_SCALAR;
  const PylithScalar fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
  };
  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const PylithScalar basis[numQuadPts*numBasis] = {
    1.0, 1.0,
    1.0, 1.0,
  };
  const PylithScalar basisDerivRef[numQuadPts*numBasis*cellDim] = {
    1.0, 1.0, 1.0, 1.0,
  };
  const PylithScalar quadPtsRef[numQuadPts*cellDim] = {
    0.0, 0.0,
  };
  const PylithScalar quadWts[numQuadPts] = { 1.5, 0.5 };
  const PylithScalar minJacobian = 1.0;

  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    topology::FieldBase::SCALAR;
  const int fiberDimE = 1;
  const PylithScalar fieldValuesE[numBasis] = {
    (1.5*1.1 + 0.5*1.2)/2.0,
    (1.5*2.1 + 0.5*2.2)/2.0,
  };
  const PylithScalar fieldScale = 2.0;

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  topology::Mesh submesh(mesh, group);

  // Set cell field
  topology::Field field(submesh);
  field.newSection(topology::FieldBase::FACES_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.label(label.c_str());
  field.scale(fieldScale);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  topology::VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    const PetscInt off = fieldVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(c));

    for(PetscInt d = 0; d < fiberDim; ++d, ++index) {
      fieldArray[off+d] = fieldValues[index];
    } // for
  } // for

  feassemble::Quadrature quadrature;
  quadrature.initialize(basis, numQuadPts, numBasis,
			basisDerivRef, numQuadPts, numBasis, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  CellFilterAvg filter;
  filter.quadrature(&quadrature);

  const topology::Field& fieldF = filter.filter(field);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldF.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldF.label()));
  CPPUNIT_ASSERT_EQUAL(fieldScale, fieldF.scale());

  topology::VecVisitorMesh fieldFVisitor(fieldF);
  const PetscScalar* fieldFArray = fieldFVisitor.localArray();CPPUNIT_ASSERT(fieldFArray);

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt c = cStart, index = 0; c < cEnd; ++c) {
    const PetscInt off = fieldFVisitor.sectionOffset(c);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldFVisitor.sectionDof(c));
    for(PetscInt d = 0; d < fiberDimE; ++d, ++index) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldFArray[off+d]/fieldValuesE[index], tolerance);
    }
  } // for

  PYLITH_METHOD_END;
} // testFilterSubMesh


// End of file 
