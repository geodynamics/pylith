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

#include "TestFieldSubMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldSubMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldSubMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 2;
      const int ncorners = 3;
      const int cells[] = {
	0, 1, 3,
	0, 3, 2,
      };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
      const char* label = "bc";
      const int groupSize = 3;
      const int groupVertices[] = {
	1, 2, 3
      };
      const int submeshVertices[] = {
	3, 4, 5
      };
    } // _TestFieldSubMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldSubMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldSubMesh::testSection(void)
{ // testSection
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(section);
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldSubMesh::testMesh(void)
{ // testMesh
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  const SubMesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldSubMesh::cellDim-1, mesh2.dimension());  
} // testMesh

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldSubMesh::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldSubMesh::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldSubMesh::testNewSection(void)
{ // testNewSection
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  field.newSection();
  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(section);
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldSubMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldSubMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldSubMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
    
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  // Create field with atlas to use to create new field
  Field<SubMesh> fieldSrc(submesh);
  fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  fieldSrc.allocate();

  const int fiberDim2 = 4;
  Field<SubMesh> field(submesh);
  field.newSection(fieldSrc, fiberDim2);
  field.allocate();

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim2, fieldVisitor.sectionDof(v));
  } // for
} // testNewSectionChart

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldSubMesh::testCloneSection(void)
{ // testCloneSection
  const int fiberDim = 3;
  const int nconstraints[] = { 0, 2, 1, 3 };
  const int constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscErrorCode err = 0;
  PetscInt vStart, vEnd;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Create field with atlas to use to create new field
  Field<SubMesh> fieldSrc(submesh);
  { // Setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    PetscSection section = fieldSrc.petscSection();CPPUNIT_ASSERT(section);
    int iV=0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
    } // for
    fieldSrc.allocate();

    int index = 0;
    iV = 0;
    for(PetscInt v = vStart; v < vEnd; ++v, index += nconstraints[iV++]) {
      err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);CHECK_PETSC_ERROR(err);
    } // for
    fieldSrc.zero();
  } // Setup source field

  Field<SubMesh> field(submesh);
  field.cloneSection(fieldSrc);
  PetscSection section = field.petscSection();CPPUNIT_ASSERT(section);
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);

  int iV = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(section, v, &cdof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], cdof);
  } // for
} // testCloneSection

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldSubMesh::testClear(void)
{ // testClear
  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);
  Field<SubMesh> field(submesh);

  field.scale(2.0);
  field.vectorFieldType(Field<SubMesh>::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), field._metadata["default"].scale);
  CPPUNIT_ASSERT_EQUAL(Field<SubMesh>::OTHER, field._metadata["default"].vectorFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._metadata["default"].dimsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldSubMesh::testAllocate(void)
{ // testAllocate
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for
} // testAllocate

// ----------------------------------------------------------------------
// Test zero().
void
pylith::topology::TestFieldSubMesh::testZero(void)
{ // testZero
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.zero();

  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } // for
  } // for
} // testZero

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldSubMesh::testComplete(void)
{ // testComplete
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> field(submesh);
  { // setup field
    field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // setup field

  field.complete();

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for
} // testComplete

// ----------------------------------------------------------------------
// Test copy().
void
pylith::topology::TestFieldSubMesh::testCopy(void)
{ // testCopy
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> fieldSrc(submesh);
  { // setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // setup source field

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.copy(fieldSrc);

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++], fieldArray[off+d], tolerance);
    } // for
  } // for
} // testCopy

// ----------------------------------------------------------------------
// Test operator+=().
void
pylith::topology::TestFieldSubMesh::testOperatorAdd(void)
{ // testOperateAdd
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesA[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };
  const PylithScalar valuesB[] = {
    10.1, 20.2, 30.3,
    10.2, 20.3, 30.4,
    10.3, 20.4, 30.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> fieldSrc(submesh);
  { // setup source field
    fieldSrc.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesA[i++];
    } // for
  } // setup source field

  Field<SubMesh> field(submesh);
  { // setup destination field
    field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesB[i++];
    } // for
  } // setup destination field

  field += fieldSrc;

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d, ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesA[i] + valuesB[i], fieldArray[off+d], tolerance);
    } // for
  } // for
} // testOperateAdd

// ----------------------------------------------------------------------
// Test dimensionalize().
void
pylith::topology::TestFieldSubMesh::testDimensionalize(void)
{ // testDimensionalize
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> field(submesh);
  { // setup field
    field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // setup field

  field.scale(scale);
  field.addDimensionOkay(true);
  field.dimensionalize();

  const PylithScalar tolerance = 1.0e-6;
  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i++]*scale, fieldArray[off+d], tolerance);
    } // for
  } // for
} // testDimensionalize

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldSubMesh::testView(void)
{ // testView
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for

  field.view("Testing view");
} // testView

// ----------------------------------------------------------------------
// Test createScatter().
void
pylith::topology::TestFieldSubMesh::testCreateScatter(void)
{ // testCreateScatter
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  const int sizeE = (vEnd-vStart) * fiberDim;  
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Make sure we can do multiple calls to createScatter().
  field.createScatter(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatter(submesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<SubMesh> field2(submesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<SubMesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<SubMesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);
} // testCreateScatter

// ----------------------------------------------------------------------
// Test createScatterWithBC().
void
pylith::topology::TestFieldSubMesh::testCreateScatterWithBC(void)
{ // testCreateScatterWithBC
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  int size = 0;
  VecGetSize(sinfo.vector, &size);
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Make sure we can do multiple calls to createScatterWithBC().
  field.createScatterWithBC(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatterWithBC(submesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<SubMesh> field2(submesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<SubMesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<SubMesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);
} // testCreateScatterWithBC

// ----------------------------------------------------------------------
// Test vector().
void
pylith::topology::TestFieldSubMesh::testVector(void)
{ // testVector
  const int fiberDim = 3;

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  Field<SubMesh> field(submesh);
  field.newSection(FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(submesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<SubMesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();
  const PetscVec vec = field.vector();
  CPPUNIT_ASSERT_EQUAL(sinfo.vector, vec);
  int size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
} // testVector

// ----------------------------------------------------------------------
// Test scatterSectionToVector().
void
pylith::topology::TestFieldSubMesh::testScatterSectionToVector(void)
{ // testScatterSectionToVector
  const char* context = "abc";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> field(submesh);
  { // setup field
    field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesE[i++];
    } // for
  } // setup field

  field.createScatter(submesh, context);
  field.scatterSectionToVector(context);

  PetscErrorCode err = 0;
  const PetscVec vec = field.vector(context);CPPUNIT_ASSERT(vec);
  PetscInt size = 0;
  err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  PetscScalar* valuesVec = NULL;
  err = VecGetArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);

  const PylithScalar tolerance = 1.0e-06;
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i], valuesVec[i], tolerance);
  err = VecRestoreArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);
} // testScatterSectionToVector

// ----------------------------------------------------------------------
// Test scatterVectorToSection().
void
pylith::topology::TestFieldSubMesh::testScatterVectorToSection(void)
{ // testScatterVectorToSection
  const char* context = "abcd";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
  };

  Mesh mesh;
  SubMesh submesh;
  _buildMesh(&mesh, &submesh);

  PetscDM dmMesh = submesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<SubMesh> field(submesh);
  field.newSection(Field<SubMesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.createScatter(submesh, context);

  const PetscVec vec = field.vector(context);CPPUNIT_ASSERT(vec);
  PetscInt size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
  PetscScalar* valuesVec = NULL;
  err = VecGetArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);
  const int sizeE = (vEnd-vStart) * fiberDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  for (int i=0; i < sizeE; ++i)
    valuesVec[i] = valuesE[i];
  err = VecRestoreArray(vec, &valuesVec);CHECK_PETSC_ERROR(err);

  field.scatterVectorToSection(context);

  const PylithScalar tolerance = 1.0e-06;
  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
    for (int iDim=0; iDim < fiberDim; ++iDim) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i++], fieldArray[off+iDim], tolerance);
    } // for
  } // for

} // testScatterVectorToSection

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldSubMesh::_buildMesh(Mesh* mesh,
					       SubMesh* submesh)
{ // _buildMesh
  assert(mesh);
  assert(submesh);

  const int cellDim = _TestFieldSubMesh::cellDim;
  const int ncells = _TestFieldSubMesh::ncells;
  const int* cells = _TestFieldSubMesh::cells;
  const int nvertices = _TestFieldSubMesh::nvertices;
  const int ncorners = _TestFieldSubMesh::ncorners;
  const int spaceDim = _TestFieldSubMesh::cellDim;
  const PylithScalar* coordinates = _TestFieldSubMesh::coordinates;
  const bool interpolate = false;

  PetscErrorCode err = 0;

  mesh->createDMMesh(_TestFieldSubMesh::cellDim);
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);  
  err = DMPlexSetChart(dmMesh, 0, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < ncells; ++c) {
    err = DMPlexSetConeSize(dmMesh, c, ncorners);CHECK_PETSC_ERROR(err);
  } // for
  err = DMSetUp(dmMesh);CHECK_PETSC_ERROR(err);
  PetscInt *cone = new PetscInt[ncorners];
  for(PetscInt c = 0; c < ncells; ++c) {
    for(PetscInt v = 0; v < ncorners; ++v) {
      cone[v] = cells[c*ncorners+v]+ncells;
    } // for
    err = DMPlexSetCone(dmMesh, c, cone);CHECK_PETSC_ERROR(err);
  } // for
  delete[] cone; cone = 0;
  err = DMPlexSymmetrize(dmMesh);CHECK_PETSC_ERROR(err);
  err = DMPlexStratify(dmMesh);CHECK_PETSC_ERROR(err);

  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  PetscScalar *coords = NULL;
  PetscInt coordSize = 0;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, ncells, ncells+nvertices);CHECK_PETSC_ERROR(err);
  for(PetscInt v = ncells; v < ncells+nvertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  } // for
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(mesh->comm(), &coordVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = 0; v < nvertices; ++v) {
    PetscInt off;
    err = PetscSectionGetOffset(coordSection, v+ncells, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = coordinates[v*spaceDim+d];
    } // for
  } // for
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(dmMesh, coordVec);CHECK_PETSC_ERROR(err);

  const int numPoints = _TestFieldSubMesh::groupSize;
  for(PetscInt i = 0; i < numPoints; ++i) {
    err = DMPlexSetLabelValue(dmMesh, _TestFieldSubMesh::label, ncells+_TestFieldSubMesh::groupVertices[i], 1);CHECK_PETSC_ERROR(err);
  } // for

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  submesh->createSubMesh(*mesh, _TestFieldSubMesh::label);
} // _buildMesh


// End of file 
