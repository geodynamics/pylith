// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestVertexFilterVecNorm.hh" // Implementation of class methods

#include "pylith/meshio/VertexFilterVecNorm.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Field.hh" // USES Field

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestVertexFilterVecNorm );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestVertexFilterVecNorm::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  VertexFilterVecNorm filter;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filter()
void
pylith::meshio::TestVertexFilterVecNorm::testFilter(void)
{ // testFilter
  PYLITH_METHOD_BEGIN;

  const char* filename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int nvertices = 4;
  const std::string label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::VECTOR;
  const PylithScalar fieldValues[nvertices*fiberDim] = {
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
  const PylithScalar fieldScale = 4.0;

  MeshIOAscii iohandler;
  topology::Mesh mesh;
  iohandler.filename(filename);
  iohandler.read(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  topology::Field field(mesh);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.vectorFieldType(fieldType);
  field.setLabel(label.c_str());
  field.scale(fieldScale);

  { // Setup vertex field
    topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    CPPUNIT_ASSERT_EQUAL(nvertices, verticesStratum.size());

    for(PetscInt v = vStart, index = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));

      for(PetscInt d = 0; d < fiberDim; ++d, ++index) {
	fieldArray[off+d] = fieldValues[index];
      } // for
    } // for
  } // Setup vertex field

  VertexFilterVecNorm filter;
  const topology::Field& fieldNorm = filter.filter(field);

  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldNorm.vectorFieldType());
  CPPUNIT_ASSERT_EQUAL(label, std::string(fieldNorm.getLabel()));
  CPPUNIT_ASSERT_EQUAL(fieldScale, fieldNorm.scale());

  topology::VecVisitorMesh fieldNormVisitor(fieldNorm);
  const PetscScalar* fieldNormArray = fieldNormVisitor.localArray();

  const PylithScalar tolerance = 1.0e-06;
  for(PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = fieldNormVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(1, fieldNormVisitor.sectionDof(v));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, fieldNormArray[off]/fieldValuesE[index++], tolerance);
  } // for

  PYLITH_METHOD_END;
} // testFilter


// End of file 
