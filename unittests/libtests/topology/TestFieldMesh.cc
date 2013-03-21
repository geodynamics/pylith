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

#include "TestFieldMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldMesh );

// ----------------------------------------------------------------------
namespace pylith {
  namespace topology {
    namespace _TestFieldMesh {
      const int cellDim = 2;
      const int nvertices = 4;
      const int ncells = 1;
      const int ncorners = 4;
      const int cells[] = { 0, 1, 2, 3 };
      const PylithScalar coordinates[] = {
	0.0, 0.0,
	1.0, 0.0,
	0.0, 1.0,
	1.0, 1.0,
      };
    } // _TestFieldMesh
  } // topology
} // pylith

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  Field<Mesh> field(mesh);
} // testConstructor

// ----------------------------------------------------------------------
// Test section().
void
pylith::topology::TestFieldMesh::testSection(void)
{ // testSection
  Mesh mesh;
  Field<Mesh> field(mesh);

  PetscSection section = field.petscSection();
  CPPUNIT_ASSERT(!section);
} // testSection

// ----------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void)
{ // testMesh
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  const Mesh& mesh2 = field.mesh();
  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, mesh2.dimension());  
} // testMesh

// ----------------------------------------------------------------------
// Test label().
void
pylith::topology::TestFieldMesh::testLabel(void)
{ // testLabel
  const std::string label = "velocity";

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  field.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(field.label()));
} // testLabel

// ----------------------------------------------------------------------
// Test vectorFieldType().
void
pylith::topology::TestFieldMesh::testVectorFieldType(void)
{ // testVectorFieldType
  const std::string label = "default";

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  field.vectorFieldType(FieldBase::SCALAR);
  CPPUNIT_ASSERT_EQUAL(FieldBase::SCALAR, field._metadata[label].vectorFieldType);
} // testVectorFieldType

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestFieldMesh::testScale(void)
{ // testScale
  const std::string label = "default";

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  const PylithScalar scale = 2.0;
  field.scale(scale);
  const PylithScalar tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, field._metadata[label].scale, tolerance);
} // testScale

// ----------------------------------------------------------------------
// Test addDimensionOkay().
void
pylith::topology::TestFieldMesh::testAddDimensionOkay(void)
{ // testAddDimensionOkay
  const std::string label = "default";

  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  CPPUNIT_ASSERT_EQUAL(false, field._metadata[label].dimsOkay);
  field.addDimensionOkay(true);
  CPPUNIT_ASSERT_EQUAL(true, field._metadata[label].dimsOkay);
} // testAddDimensionOkay

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestFieldMesh::testSpaceDim(void)
{ // testSpaceDim
  Mesh mesh;
  _buildMesh(&mesh);
  Field<Mesh> field(mesh);

  CPPUNIT_ASSERT_EQUAL(_TestFieldMesh::cellDim, field.spaceDim());
} // testSpaceDim

// ----------------------------------------------------------------------
// Test newSection().
void
pylith::topology::TestFieldMesh::testNewSection(void)
{ // testNewSection
  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  field.newSection();
  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscSection section = field.petscSection();CPPUNIT_ASSERT(section);
} // testNewSection

// ----------------------------------------------------------------------
// Test newSection(points).
void
pylith::topology::TestFieldMesh::testNewSectionPoints(void)
{ // testNewSectionPoints
  const int fiberDim = 2;
  const std::string& label = "field A";

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  field.label(label.c_str());
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionPoints

// ----------------------------------------------------------------------
// Test newSection(int_array).
void
pylith::topology::TestFieldMesh::testNewSectionPointsArray(void)
{ // testNewSectionPointsArray
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const int npts = (vEnd-vStart) / 2;
  int_array pointsIn(npts);
  int_array pointsOut(vEnd-vStart - npts);
  int count = 0;
  size_t iIn = 0;
  size_t iOut = 0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    if (count % 2  == 0)
      pointsIn[iIn++] = v;
    else
      pointsOut[iOut++] = v;
    ++count;
  } // for
  CPPUNIT_ASSERT_EQUAL(iIn, pointsIn.size());
  CPPUNIT_ASSERT_EQUAL(iOut, pointsOut.size());

  field.newSection(pointsIn, fiberDim);
  field.allocate();

  // Points in array should have a fiber dimension of fiberDim.
  VecVisitorMesh fieldVisitor(field);
  for(int i = 0; i < pointsIn.size(); ++i) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(pointsIn[i]));
  } // for

  // Points not int array should have a fiber dimension of zero.
  PetscInt pStart, pEnd;
  PetscSection section = field.petscSection();CPPUNIT_ASSERT(section);
  PetscErrorCode err = PetscSectionGetChart(section, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  for(int i = 0; i < pointsOut.size(); ++i) {
    if (pointsOut[i] >= pStart && pointsOut[i] < pEnd) {
      CPPUNIT_ASSERT_EQUAL(0, fieldVisitor.sectionDof(pointsOut[i]));
    } // if
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionPointsArray

// ----------------------------------------------------------------------
// Test newSection(domain).
void
pylith::topology::TestFieldMesh::testNewSectionDomain(void)
{ // testNewSectionDomain
  const int fiberDim = 2;

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionDomain

// ----------------------------------------------------------------------
// Test newSection(field).
void
pylith::topology::TestFieldMesh::testNewSectionField(void)
{ // testNewSectionField
  const int fiberDim = 3;
    
  Mesh mesh;
  _buildMesh(&mesh);

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  fieldSrc.allocate();

  const int fiberDim2 = 5;
  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(fieldSrc, fiberDim2);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  VecVisitorMesh fieldVisitor(field);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    CPPUNIT_ASSERT_EQUAL(fiberDim2, fieldVisitor.sectionDof(v));
  } // for

  const char *name = NULL;
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  PetscErrorCode err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testNewSectionField

// ----------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCloneSection(void)
{ // testCloneSection
  const PetscInt fiberDim = 3;
  const PetscInt nconstraints[] = { 0, 2, 1, 3 };
  const PetscInt constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    PetscSection section = fieldSrc.petscSection();CPPUNIT_ASSERT(section);
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
    } // for
    fieldSrc.allocate();

    int index = 0;
    for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
      err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);CHECK_PETSC_ERROR(err);
    }
    fieldSrc.zero();
    fieldSrc.createScatter(mesh);
    fieldSrc.createScatter(mesh, "A");
  } // Setup source field

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.cloneSection(fieldSrc);
  PetscSection section = field.petscSection();
  PetscVec vec     = field.localVector();
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  for(PetscInt v = vStart, iV = 0; v < vEnd; ++v) {
    PetscInt dof, cdof;
    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetConstraintDof(section, v, &cdof);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDim, dof);
    CPPUNIT_ASSERT_EQUAL(nconstraints[iV++], cdof);
  } // for

  // Verify vector scatters were also copied.
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters[""].dm,  field._scatters[""].dm);
  CPPUNIT_ASSERT_EQUAL(fieldSrc._scatters["A"].dm, field._scatters["A"].dm);
  const char *name = NULL;
  err = PetscObjectGetName((PetscObject) vec, &name);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(label, std::string(name));
} // testCloneSection

// ----------------------------------------------------------------------
// Test clear().
void
pylith::topology::TestFieldMesh::testClear(void)
{ // testClear
  Mesh mesh(_TestFieldMesh::cellDim);
  Field<Mesh> field(mesh);

  field.scale(2.0);
  field.vectorFieldType(Field<Mesh>::TENSOR);
  field.addDimensionOkay(true);
  
  field.clear();

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), field._metadata["default"].scale);
  CPPUNIT_ASSERT_EQUAL(Field<Mesh>::OTHER, field._metadata["default"].vectorFieldType);
  CPPUNIT_ASSERT_EQUAL(false, field._metadata["default"].dimsOkay);
} // testClear

// ----------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void)
{ // testAllocate
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
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
pylith::topology::TestFieldMesh::testZero(void)
{ // testZero
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
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
// Test zero().
void
pylith::topology::TestFieldMesh::testZeroAll(void)
{ // testZeroAll
  const PetscInt fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const PetscInt nconstraints[] = { 0, 2, 1, 3 };
  const PetscInt constraints[] = {
              // 0
    0, 2,     // 1
    2,        // 2
    0, 1, 2,  // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  // Create field and set constraint sizes
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  PetscSection section = field.petscSection();CPPUNIT_ASSERT(section);
  PetscErrorCode err = 0;
  PetscInt index = 0;
  for(PetscInt v = vStart, iV=0; v < vEnd; ++v) {
    err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
  } // for
  field.allocate();
  PetscVec vec = field.localVector();CPPUNIT_ASSERT(vec);
  for(PetscInt v = vStart, iV = 0; v < vEnd; ++v, index += nconstraints[iV++]) {
    err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[index]);CHECK_PETSC_ERROR(err);
  } // for
  field.zero();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.zeroAll();
  
  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, fieldArray[off+d], tolerance);
    } // for
  } // for
} // testZeroAll

// ----------------------------------------------------------------------
// Test complete().
void
pylith::topology::TestFieldMesh::testComplete(void)
{ // testComplete
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  VecVisitorMesh fieldVisitor(field);
  PetscScalar* fieldArray = fieldVisitor.localArray();

  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d)
      fieldArray[off+d] = valuesNondim[i++];
  } // for
  fieldVisitor.clear();

  field.complete();

  // Expect no change for this serial test
  fieldVisitor.initialize(field);
  fieldArray = fieldVisitor.localArray();
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
pylith::topology::TestFieldMesh::testCopy(void)
{ // testCopy
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesNondim[i++];
    } // for
  } // Setup source field
    
  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
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
pylith::topology::TestFieldMesh::testOperatorAdd(void)
{ // testOperateAdd
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesA[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };
  const PylithScalar valuesB[] = {
    10.1, 20.2, 30.3,
    10.2, 20.3, 30.4,
    10.3, 20.4, 30.5,
    10.4, 20.5, 30.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    fieldSrc.allocate();
    VecVisitorMesh fieldVisitor(fieldSrc);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesA[i++];
    } // for
  } // Setup source field
    

  Field<Mesh> field(mesh);
  { // Setup destination field
    field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesB[i++];
    } // for
  } // Setup destination field

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
pylith::topology::TestFieldMesh::testDimensionalize(void)
{ // testDimensionalize
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> field(mesh);
  { // setup field
    field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
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

  VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();
  const PylithScalar tolerance = 1.0e-6;
  for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
    const PetscInt off = fieldVisitor.sectionOffset(v);
    for(PetscInt d = 0; d < fiberDim; ++d, ++i) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesNondim[i]*scale, fieldArray[off+d], tolerance);
    } // for
  } // for
} // testDimensionalize

// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldMesh::testView(void)
{ // testView
  const int fiberDim = 3;
  const PylithScalar scale = 2.0;
  const PylithScalar valuesNondim[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
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
pylith::topology::TestFieldMesh::testCreateScatter(void)
{ // testCreateScatter
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  const int sizeE = (vEnd-vStart) * fiberDim;
  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Check vector name
  const char* vecname = 0;
  PetscObjectGetName((PetscObject)sinfo.vector, &vecname);
  CPPUNIT_ASSERT_EQUAL(label, std::string(vecname));

  // Make sure we can do multiple calls to createScatter().
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatter(mesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);

} // testCreateScatter

// ----------------------------------------------------------------------
// Test createScatterWithBC().
void
pylith::topology::TestFieldMesh::testCreateScatterWithBC(void)
{ // testCreateScatterWithBC
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  const std::string& label = "field A";
  field.label(label.c_str());
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();

  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  const int sizeE = (vEnd-vStart) * fiberDim;
  int size = 0;
  VecGetSize(sinfo.vector, &size);
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  // Check vector name
  const char* vecname = 0;
  PetscObjectGetName((PetscObject)sinfo.vector, &vecname);
  CPPUNIT_ASSERT_EQUAL(label, std::string(vecname));

  // Make sure we can do multiple calls to createScatterWithBC().
  field.createScatterWithBC(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());

  // Create another scatter.
  field.createScatterWithBC(mesh, "B");
  CPPUNIT_ASSERT_EQUAL(size_t(2), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfoB = field._getScatter("B");
  CPPUNIT_ASSERT(sinfoB.dm);
  CPPUNIT_ASSERT(sinfoB.vector);

  Field<Mesh> field2(mesh);
  field2.cloneSection(field);
  CPPUNIT_ASSERT_EQUAL(size_t(2), field2._scatters.size());

  const Field<Mesh>::ScatterInfo& sinfo2 = field2._getScatter("");
  CPPUNIT_ASSERT(sinfo2.dm);
  CPPUNIT_ASSERT(sinfo2.vector);

  const Field<Mesh>::ScatterInfo& sinfo2B = field2._getScatter("B");
  CPPUNIT_ASSERT(sinfo2B.dm);
  CPPUNIT_ASSERT(sinfo2B.vector);

} // testCreateScatterWithBC

// ----------------------------------------------------------------------
// Test vector().
void
pylith::topology::TestFieldMesh::testVector(void)
{ // testVector
  const int fiberDim = 3;

  Mesh mesh;
  _buildMesh(&mesh);

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  
  CPPUNIT_ASSERT_EQUAL(size_t(0), field._scatters.size());
  field.createScatter(mesh);
  CPPUNIT_ASSERT_EQUAL(size_t(1), field._scatters.size());
  const Field<Mesh>::ScatterInfo& sinfo = field._getScatter("");
  CPPUNIT_ASSERT(sinfo.dm);
  CPPUNIT_ASSERT(sinfo.vector);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
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
pylith::topology::TestFieldMesh::testScatterSectionToVector(void)
{ // testScatterSectionToVector
  const char* context = "abc";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> field(mesh);
  { // setup field
    field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
    field.allocate();
    VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();
    for(PetscInt v = vStart, i = 0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      for(PetscInt d = 0; d < fiberDim; ++d)
	fieldArray[off+d] = valuesE[i++];
    } // for
  } // setup field

  field.createScatter(mesh, context);
  field.scatterSectionToVector(context);

  const PetscVec vec = field.vector(context);CPPUNIT_ASSERT(vec);
  PetscInt size = 0;
  PetscErrorCode err = VecGetSize(vec, &size);CHECK_PETSC_ERROR(err);
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
pylith::topology::TestFieldMesh::testScatterVectorToSection(void)
{ // testScatterVectorToSection
  const char* context = "abcd";
  const int fiberDim = 3;
  const PylithScalar valuesE[] = {
    1.1, 2.2, 3.3,
    1.2, 2.3, 3.4,
    1.3, 2.4, 3.5,
    1.4, 2.5, 3.6,
  };

  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  Field<Mesh> field(mesh);
  field.newSection(Field<Mesh>::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.createScatter(mesh, context);

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
// Test splitDefault().
void
pylith::topology::TestFieldMesh::testSplitDefault(void)
{ // testSplitDefault
  const PetscInt spaceDim = _TestFieldMesh::cellDim;
  const PetscInt numFields = spaceDim;
  const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
  const PetscInt constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with section to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    for(PetscInt f = 0; f < numFields; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.addField(msg.str().c_str(), 1);
    } // for
    fieldSrc.setupFields();
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    for(PetscInt f = 0; f < spaceDim; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.updateDof(msg.str().c_str(), Field<Mesh>::VERTICES_FIELD, 1);
    } // for
    PetscSection section = fieldSrc.petscSection();CPPUNIT_ASSERT(section);
    PetscInt iV = 0, iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV];
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int field = constraints[iC++];
        err = PetscSectionAddFieldConstraintDof(section, v, field, 1);CHECK_PETSC_ERROR(err);
      } // for
    } // for
    fieldSrc.allocate();

    iV = 0; iC = 0;
    PetscInt zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV++];
      if (nconstraintsVertex > 0) {
        err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]);CHECK_PETSC_ERROR(err);
      } // if
      for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const PetscInt field = constraints[iC++];
        err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero);CHECK_PETSC_ERROR(err);
      } // for
    } // for
  } // Setup source field

  PetscSection section = fieldSrc.petscSection();CPPUNIT_ASSERT(section);
  PetscInt numSectionFields;
  err = PetscSectionGetNumFields(section, &numSectionFields);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

  for(PetscInt f = 0; f < numFields; ++f) {
    PetscInt iV = 0, iC = 0;

    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      PetscInt fdof, fcdof;
      PetscBool isConstrained = PETSC_FALSE;
      const PetscInt nconstraintsVertex = nconstraints[iV];

      err = PetscSectionGetFieldDof(section, v, f, &fdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(1, fdof);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == f) isConstrained = PETSC_TRUE;
      const PetscInt constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE, fcdof);
    } // for
  } // for
} // testSplitDefault

// ----------------------------------------------------------------------
// Test cloneSection() with split field.
void
pylith::topology::TestFieldMesh::testCloneSectionSplit(void)
{ // testCloneSectionSplit
  const PetscInt spaceDim = _TestFieldMesh::cellDim;
  const PetscInt numFields = spaceDim;
  const PetscInt nconstraints[4] = { 1, 2, 0, 1 };
  const PetscInt constraints[4] = {
    1,     // 0
    0, 1,  // 1
           // 2
    0,     // 3
  };
    
  Mesh mesh;
  _buildMesh(&mesh);

  PetscDM dmMesh = mesh.dmMesh();CPPUNIT_ASSERT(dmMesh);
  Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();

  PetscErrorCode err = 0;

  // Create field with atlas to use to create new field
  Field<Mesh> fieldSrc(mesh);
  { // Setup source field
    for(PetscInt f = 0; f < numFields; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.addField(msg.str().c_str(), 1);
    } // for
    fieldSrc.setupFields();
    fieldSrc.newSection(Field<Mesh>::VERTICES_FIELD, spaceDim);
    for(PetscInt f = 0; f < spaceDim; ++f) {
      std::ostringstream msg;
      msg << "Field "<<f;
      fieldSrc.updateDof(msg.str().c_str(), Field<Mesh>::VERTICES_FIELD, 1);
    } // for
    PetscSection section = fieldSrc.petscSection();
    CPPUNIT_ASSERT(section);
    PetscInt iV = 0, iC = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV];
      err = PetscSectionAddConstraintDof(section, v, nconstraints[iV++]);CHECK_PETSC_ERROR(err);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const int field = constraints[iC++];
        err = PetscSectionAddFieldConstraintDof(section, v, field, 1);CHECK_PETSC_ERROR(err);
      } // for
    } // for
    fieldSrc.allocate();

    iV = 0; iC = 0;
    PetscInt zero = 0;
    for(PetscInt v = vStart; v < vEnd; ++v) {
      const int nconstraintsVertex = nconstraints[iV++];
      if (nconstraintsVertex > 0) {
        err = PetscSectionSetConstraintIndices(section, v, (PetscInt *) &constraints[iC]);CHECK_PETSC_ERROR(err);
      } // if
      for(PetscInt iConstraint=0; iConstraint < nconstraintsVertex; ++iConstraint) {
        const PetscInt field = constraints[iC++];
        err = PetscSectionSetFieldConstraintIndices(section, v, field, &zero);CHECK_PETSC_ERROR(err);
      } // for
    } // for
  } // Setup source field

  Field<Mesh> field(mesh);
  field.cloneSection(fieldSrc);

  PetscSection section = field.petscSection();CPPUNIT_ASSERT(section);
  PetscInt numSectionFields;
  err = PetscSectionGetNumFields(section, &numSectionFields);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(numFields, numSectionFields);

  for(PetscInt f = 0; f < numFields; ++f) {
    PetscInt iV = 0, iC = 0;

    for(PetscInt v = vStart; v < vEnd; ++v, ++iV) {
      PetscInt fdof, fcdof;
      PetscBool isConstrained = PETSC_FALSE;
      const PetscInt nconstraintsVertex = nconstraints[iV];

      err = PetscSectionGetFieldDof(section, v, f, &fdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetFieldConstraintDof(section, v, f, &fcdof);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(1, fdof);
      for(PetscInt iConstraint = 0; iConstraint < nconstraintsVertex; ++iConstraint)
        if (constraints[iC++] == f) isConstrained = PETSC_TRUE;
      const PetscInt constraintDimE = (!isConstrained) ? 0 : 1;
      CPPUNIT_ASSERT_EQUAL(constraintDimE, fcdof);
    } // for
  } // for
} // testCloneSectionSplit

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_buildMesh(Mesh* mesh)
{ // _buildMesh
  assert(mesh);

  const int cellDim = _TestFieldMesh::cellDim;
  const int ncells = _TestFieldMesh::ncells;
  const int* cells = _TestFieldMesh::cells;
  const int nvertices = _TestFieldMesh::nvertices;
  const int ncorners = _TestFieldMesh::ncorners;
  const int spaceDim = _TestFieldMesh::cellDim;
  const PylithScalar* coordinates = _TestFieldMesh::coordinates;
  const bool interpolate = false;

  PetscErrorCode err = 0;

  mesh->createDMMesh(_TestFieldMesh::cellDim);
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
  PetscInt coordSize;

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

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
} // _buildMesh


// End of file 
